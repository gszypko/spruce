# mhdtoy
A two-dimensional MHD simulation with support for modular physics components.

# Compiling
This code was only developed using the GNU g++ compiler. Support for other C++ compilers is not guaranteed, but the code was developed across Windows, Linux, and macOS systems and so should be usable from any of the above.

The included Makefile can be used to compile by invoking `make` from the command line. The Makefile also allows you to specify the name of the executable file with the variable `EXEC`. For instance, `make EXEC=test_run` will compile into an executable named `test_run`. The default executable name is `run`, which will be used as the name of the executable in all examples below.

# Running
`mhdtoy` can be run in three different modes, depending on the method used to specify the initial state of the simulation run.

The resulting simulation will be output as a time series to the file `mhd.out` in the run's output directory. Other files written to the output directory are:
- `init.state`: contains the full initial state of the system
- `end.state`: contains the full final state of the system
- A copy of the `.config` file used for the run
- `plasma.settings`: (if run in Problem Generator Mode) encodes the particular `settings` used for the run. Note that this is not the same as the original `.settings` file if multiple sets of `settings` are possible from that file, depending on the specified `job_index`

### Custom Input Mode

A user-generated set of grids, in the form of a `.state` file, can be specified to initialize the simulation. See below for more details on formatting a `.state` file. Running in this mode requires the following syntax:

`./run  [-g grid_file] [-c config_file] [-o out_directory] [-d time_duration]`
- `-g` or `--grids` (Required) specifies the `.state` file to be used to initialize the simulation run
- `-c` or `--config` (Required) specifies the `.config` file to be used for configuring the simulation run. 
- `-o` or `--output` (Required) specifies the name of the directory to which all output will be written. If the directory does not exist it will be created, and existing simulation outputs in that directory will be overwritten. 
- `-d` or `--duration` (Required) specifies how long the simulation is to be continued, in units of simulation time.

### Continue Mode

Any previous run can be continued from where it left off, appending to the existing output files. Running in this mode requires the following syntax:

`./run [-p prev_out_directory] [-d time_duration]`
- `-p` or `--prev` (Required) specifies the `out_directory` of the previous simulation run to be continued. Note that if the previous run had a `job_index` specified, the corresponding subdirectory must also be specified.
- `-d` or `--duration` (Required) specifies how long the simulation is to be continued, in units of simulation time.

### Problem Generator Mode (Deprecated)

The compiled executable includes a number of automated initial state generators, the specifics of which are dictated by a human-readable `.settings` file. In this mode, the initial state will be generated at run-time and the simulation will immediately proceed from there. Running in this mode requires the following syntax:

`./run [-t solar|ucnp] [-s settings_file] [-i job_index] [-c config_file] [-o out_directory] [-d time_duration]`
- `-t` or `--type` (Required) specifies the problem generator to use. Currently the only valid `type`s are `ucnp` (ultracold neutral plasma) and `solar` (solar coronal plasma). 
- `-s` or `--settings` (Required) specifies the `.settings` file containing `settings` that are fed into the problem generator. The `settings` to be specified are dictated by the `type` of the problem generator being used.
- `-i` or `--index` (Optional) specifies the index (starting at zero) of the combination of `settings` to use. The `.settings` file can be written to allow multiple values for any number of different `settings` (see below). Specifying the `index` of the run uniquely specifies a combination of those `settings` to use. Defaults to `0` if not specified. If specified, output for the run will be written to a subdirectory of `out_directory`: `out_directory/array0/`, `out_directory/array1/`, etc.
- `-c` or `--config` (Required) specifies the `.config` file to be used for configuring the simulation run. 
- `-o` or `--output` (Required) specifies the name of the directory to which all output will be written. If the directory does not exist it will be created, and existing simulation outputs in that directory will be overwritten. 
- `-d` or `--duration` (Required) specifies how long the simulation is to be continued, in units of simulation time.


# The State File

The variables that need to be specified throughout the simulation domain at run-time are `{d_x,d_y,pos_x,pos_y,rho,temp,mom_x,mom_y,be_x,be_y,bi_,bi_y,grav_x,grav_y}`. These are
- `d_x` and `d_y`: The width of each cell, in cm, in the x- and y-directions respectively. In the current implementation these need not be uniform, but they can only change along the corresponding direction. That is, `d_x` must be the same for all y-values given a corresponding x-value, and likewise for `d_y`.
- `pos_x` and `pos_y`: The positions of each cell's center, in cm, in the x- and y-directions. These must correspond to the given `d_x` and `d_y`, but the position of the origin is up to the user to decide.
- `rho`: The mass density, in g cm^-3.
- `temp`: The temperature, in K
- `mom_x` and `mom_y`: The momentum density, in erg cm^-3
- `be_x` and `be_y`: The external magnetic field, in gauss. This magnetic field is held constant, and can act on the plasma but cannot be acted on by it.
- `bi_x` and `bi_y`: The induced magnetic field, in gauss. This is the magnetic field generated by motion of the plasma. A self-consistent MHD equilibrium, for instance, should have its field specified using `bi_x` and `bi_y` while `be_x` and `be_y` are set to zero.
- `grav_x` and `grav_y`: The acceleration due to gravity at each position, in cm s^-2.

The state file must be formatted as follows:
```
xdim,ydim
[xdim],[ydim]
ion_mass
[ion mass]
adiabatic_index
[adiabatic index]
t=[start time]
[variable name]
[variable grid]
[variable name]
[variable grid]
...
```
where the bracketed text should be repalced as follows:
- `xdim`: The number of grid cells in the x-direction.
- `ydim`: The number of grid cells in the y-direction.
- `ion mass`: The ion mass to be used in the ideal MHD equations, in grams. This is used to convert between `rho` and number density `n`.
- `adiabatic index`: The adiabatic index to be used in the simulation. This is used to convert between pressure `press` and thermal energy `thermal_energy`.
- `start time`: The time, in seconds, at which to start the simulation run. Unless continuing from a previous run, this should probably be set to zero.
- `variable name`: One of the variable names given in the list above. These can be specified in any order, and any variables not specified will be set to zero.
- `variable grid`: The values of the immediately previous `variable name` specified at every position in the simulation grid. This should be formatted as a set of comma-separated values, where each row corresponds to an x-value, running from top to bottom. Each column corresponds to a y-value, running from left to right. This means that each `variable grid` should consist of `xdim` rows and `ydim` columns, with the conventional matrix indexing `(i,j)` corresponding to the x- and y-indices, respectively.

See the `example.state` file for an example of a correctly formatted `.state` file.

# The Config File

The `.config` file determines the base runtime behavior of the simulation. Each line of the `.config` file must be formatted as `config_name = config_value`. Whitespace (except for newlines) is ignored, and all text following `#` is treated as a comment. The following `config`s are applicable to all simulation runs:

### General

- `time_integrator`: One of `{euler,rk2,rk4}`. Specifies the time integration scheme used for the main time integration loop.
- `max_iterations`: Integer. Specifies the maximum number of iterations allowed in the base MHD simulation. If negative, no maximum number of iterations will be enforced.
- `iter_output_interval`: Integer. Specifies how frequently to write the state of the system to the `mhd.out` file in terms of iterations. If negative, number of iterations is ignored for the purposes of writing to `mhd.out`.
- `time_output_interval`: Decimal. Specifies how frequently to write the state of the system to the `mhd.out` file in terms of simulation time. If negative, simulation time is ignored for the purposes of writing to `mhd.out`.
- `std_out_interval`: Integer. Specifies how frequently, in terms of iterations, to write a short (one-line) update on the simulation run to standard output.
- `safe_state_mode`: One of `{true,false}`. When `true`, a full `.state` file will be written periodically during the simulation run so that it can be continued even if execution is halted unceremoniously. Subsequent state files will be written and overwritten to `mhd0.state` and `mhd1.state` in the output directory; upon successful completion they will be removed in favor of the final `mhd.state`. When `false`, a full `.state` file will only be written out upon completion of the simulation run.
- `safe_state_interval`: Integer. When `safe_state_mode` is `true`, specifies the number of iterations between writing successive `.state` files.
- `output_flags`: A comma-separated list composed of any of `{d_x,d_y,pos_x,pos_y,rho,temp,mom_x,mom_y,be_x,be_y,grav_x,grav_y,press,thermal_energy,kinetic_energy,div_bi,bi_x,bi_y,dt,v_x,v_y,n,div_be,b_hat_x,b_hat_y,b_magnitude,v_a}` and/or any outputs specified by active `module`s. Specifies the variables that are written to the `mhd.out` file whenever required.

### Boundary conditions
- `x_bound_1`, `x_bound_2`, `y_bound_1`, and `y_bound_2`: One of `{periodic, open, fixed, reflect}` for each. Specifies the lower (1) and upper (2) boundary conditions applied in the `x` and `y` directions.
- `open_boundary_strength`: Decimal. When any boundary conditions is set to `open`, controls how strongly outward advection is imposed at the boundary as a multiple of the local sound speed.
- `open_boundary_decay_base`: Decimal between 0.0 and 1.0. When any boundary conditions is set to `open`, controls the imposed outward decay of mass density and thermal energy at `open` boundaries.

### Safety factors
- `epsilon`: Decimal between 0.0 and 1.0. Specifies the epsilon used when computing the time step for the current MHD iteration. Multiplied into the upper bound given by the CFL condition to give the time step.
- `epsilon_viscous`: Decimal between 0.0 and 1.0. Specifies the strength of the artificial viscosity in the simulation.

### Variable lower bounds
- `rho_min`: Decimal. Minimum allowed value for the mass density.
- `temp_min`: Decimal. Minimum allowed value for the temperature.
- `thermal_energy_min`: Decimal. Minimum allowed value for the thermal energy.

## Module Configs

The `.config` file also allows the user to configure the behavior of any additional `modules` on top of the base MHD simulation. The proper syntax for specifying a module's `config`s is

```
module_name = true|false
{
    module_config_1 = ...
    module_config_2 = ...
    ...
}
```

which can be included anywhere in the `.config` file. When the `module` is set to `true`, it will be activated at run time and fed the specified module configs in the `{}`. When the `module` is set to `false`, it will not be activated and the module configs following it in the `{}` will be ignored. 

# Modules

`mhdtoy` is designed to allow for the implementation of `Module`s that encompass physics phenomena outside of the base MHD physics of the simulation. These Modules can be easily activated, deactivated, and reconfigured without modifying the base MHD functionality. This is meant to be the primary method for the user to customize the run-time behavior of the code.

## Basic Modular Structure

The `PlasmaDomain` object that handles the main simulation loop of the plasma maintains a `ModuleHandler` object, which is responsible for administrating the behavior of all individual `Module`s in effect. All `Module`s maintain a reference to the grandparent `PlasmaDomain` and have private `friend` access, meaning they can access and modify any members of the grandparent `PlasmaDomain` by accessing their `m_pd` object during run time.

All Modules should be implemented as a `derived` class of the abstract `Module` class, declared in `source/modules/module.hpp`. The following functions may be overridden from the base Module class; note that, unless marked as **REQUIRED**, the functions below do nothing if they aren't overridden. The **REQUIRED** functions should be overridden.

- `void iterateModule(double dt) override;`
    - Evolve the system in time by dt according to the Module. Note that calling `m_pd.propagateChanges()` is likely desired at the end of this function, and potentially between iterations of a loop within this function; this will ensure that all interconnected variables (e.g. `n` and `rho`) remain consistent.
- `void preIterateModule(double dt) override;`
    - Perform any processing that occurs before the main iteration step. All Modules are guaranteed to run this function before any of them run `iterateModule()` or `postIterateModule()`.
- `void postIterateModule(double dt) override;`
    - Perform any processing that occurs after the main iteration step. All Modules are guaranteed to run this function after all of them run `iterateModule()` and `preIterateModule()`.
- `std::string commandLineMessage() const override;`
    - Returns a short message to print to stdout related to the most recent iteration of the Module. Should not include any line breaks.
- `void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const override;`
    - Returns data to write to the `PlasmaDomain`'s `mhd.out` file corresponding to the most recent iteration of the Module. Any variables to output must have their names appended to `var_names` and the corresponding Grids appended to `var_grids`, in the same order. 
- `void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;` (**REQUIRED**)
    - Applies the module configs from the current run's `.config` file. `lhs` contains the names of the module configs for the current Module from the `.config` file, and `rhs` contains the corresponding value(s) in the same order, as strings.
