# mhdtoy
A two-dimensional quasi-MHD simulation for static magnetic field

# Compiling
The Makefile is very dumb and doesn't keep track of any dependencies. It only has phony targets clang, intel, and gnu corresponding to full compilation commands using clang++, icc, and g++ compilers respectively. For instance, execute `make intel` to compile using icc.

# Running
Usage: `./mhdtoy [-o run_name] [-s settings_file] [-S state_file] [-i job_index]`
- `-o` or `--output` allows you to specify the name of the run, which is used as the prefix for the resulting `.out` and `.state` files. Default is `output`.
- `-s` or `--settings` allows you to specify the settings file to be used for configuring the simulation. Default is `default.settings`.
- `-S` or `--state` allows you to specify a .state file from a previous run to be used for initialization. If none specified, a default initialization is used.
- `-i` or `--index` allows you to specify the index of job array index of the current run (for use with Slurm job arrays, see below). Default value is 0.

The resulting simulation will be output as a time series to the file `[run_name].out` and as a state file to `[run_name].state` (which can be used to initialize a future simulation from the end point of that simulation).

# Simulation Settings
All settings are set at run time by the contents of the specified `.settings` file. See `default.settings` for the required syntax of the `.settings` file.

For all settings except for `output_flags` and `state_flags`, specifying multiple comma-separated settings in a single line allows multiple settings configurations to be indexed automatically with a job array. The appropriate array indices range from zero to the total number of settings combinations possible from the `.settings` file minus one. If not specified, the job array index defaults to zero, which will take the first setting on each line. The `run_name` is appended to indicate the settings chosen, if multiple combinations are possible from the `.settings` file.

For example, suppose the `.settings` file contains two lines with multiple values specified: `b_0 = 100, 10, 1` and `radiative_losses = true, false`. This file has six possible combinations of values for `b_0` and `radiative_losses`, meaning `job_index` can range from 0 to 5. If run with `-i 4`, for instance, the simulation will run with `b_0 = 10` and `radiative_losses = false` (along will all of the other single-value settings specified in the file), and `run_name` will be appended to `[run_name]-b_0:10-radiative_losses:false`. (Note: the ordering of setting combinations by index is dependent on the order of the multi-valued setting lines in the `.settings` file.)

Note that specifying multiple values for `output_flags` or `state_flags` will not be treated as multiple possible setting combinations; instead, all of the given variables will be given `true` flags for `.out` and `.state` files, respectively.

Also, note that the `xdim`, `ydim`, `dx`, and `dy` in a chosen `.state` file will override the values specified in the `.settings` file.
