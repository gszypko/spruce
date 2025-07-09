
# mhs_terradas

A set of scripts for generating 2D magnetohydrostatic (MHS) equilibria, using the technique described by Terradas et al. (2022) https://www.aanda.org/10.1051/0004-6361/202142975 . This is split between three different Python scripts, all written to be executed from the command line.

## The Scripts

`compute.py` does the heavy duty calculation of a normalized version of the MHS equilibrium. See the script for details of this normalization. You can modify the constants at the top of the script to control things like domain size and dimensions, resolution, locations and relative magnitudes and signs of magnetic poles, and the size parameters for the lower-boundary profile. This creates a normalized.txt file.
```
Usage: ./compute.py output_directory
output_directory is a folder to put all output into (normalized.txt and a few plots),
which will be created if it doesn't already exist
```

`denormalize.py` takes the normalized.txt file and applies the desired physical parameters to create a "de-normalized" version of the MHS equilibrium, in physical units. You can modify the constants at the top of the script to control things like the reference plasma beta, coronal temperature, and absolute magnetic field strength/sign. This produces an `.state` file.

```
Usage: ./denormalize.py input_folder output_folder [output_filename]
input_folder is the path to a folder containing a "normalized.txt" file representing a solution from mhs_terradas/compute.py
output_folder is the path to a folder to place the final .state file (creates the folder if it doesn't exist)
output_filename, optionally, is the name of the output file that will be saved to output_folder; defaults to mhs.state
```

`denormalize_join.py` does the same as `denormalize.py`, but instead takes two different normalized.txt files for a closed- and open-field solution and superimposes the two. We use this to generate a scenario at the boundary of an open-field and closed-field feature.

```
Usage: ./denormalize_join.py AR_directory CH_directory output_folder [output_filename]
AR_folder is the path to a folder containing a "normalized.txt" file representing a closed-field solution
CH_folder is the path to a folder containing a "normalized.txt" file representing an open-field solution
output_folder is the path to a folder to place the final .state file (creates the folder if it doesn't exist)
output_filename, optionally, is the name of the output file that will be saved to output_folder; defaults to mhs.state
```

  ## Output File Format

The final outputs from `denormalize.py` or `denormalize_join.py` are basically just CSV files (named `mhs.state` by default), with all variables stored with row-major indexing (that is, x-value changes with row and y-value changes with column) following a single line containing the name of the variable. They also contain several header comment lines, delineated by `#`. It should hopefully be straightforward to convert this output to your desired file format, or otherwise to modify the scripts to output directly to your desired format. The variables written out are mass density "rho", temperature "temp", momentum density "mom_x", "mom_y", "mom_z", curl-free magnetic field "be_x","be_y","be_z", non-curl-free magnetic field "bi_x","bi_y","bi_z", gravitational acceleration "grav_x","grav_y", x- and y-positions of the grid cell centers "pos_x","pos_y", and x- and y size of the grid cells "d_x","d_y". The total magnetic field is be+bi, but I leave them separate to make some calculations in my simulation simpler. All denormalized variables are in Gaussian CGS units.

  ## Limitations

A few key limitations to this method, in the interest of completeness:

- This is strictly a 2D method, and I don't believe there's any straightforward way to extend the approach into generating 3D equilibria. When generating these states I also add a small, homogeneous, nonzero magnetic field component in the out-of-plane direction (z) in order to allow the system to evolve in 2.5D.

- This method works only in the low-beta regime; the authors make clear that the solutions generated should be considered properly coronal and not chromospheric. The core approximation is a first-order Taylor expansion that applies a small perturbation to a potential magnetic field that counteracts pressure gradients in the plasma. The expansion parameter is proportional to the plasma beta, therefore this first-order perturbation is only a good approximation for low-beta plasmas (when the perturbation away from a potential solution is small). It is possible, as the authors mention, to augment this approach to a second-, third-, or higher-order approximation instead (see Eq. 39 and surrounding discussion) which should grant more leniency in the plasma beta, but I have not looked into this in detail or implemented it.

- In the current implementation gravity is constant everywhere, equal to the gravitational acceleration at the solar surface.

- A minor detail, but the paper and my scripts refer to "coronal holes" and "active regions" but the technique is rather flexible and can easily create states that are open- or closed-field but wouldn't really qualify as a coronal hole or active region-- "CH" and "AR" make for a useful shorthand.