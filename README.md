# mhdtoy
A two-dimensional quasi-MHD simulation for static magnetic field

# Compiling
The Makefile is very dumb and doesn't keep track of any dependencies. It only has phony targets clang, intel, and gnu corresponding to full compilation commands using clang++, icc, and g++ compilers respectively. For instance, execute `make intel` to compile using icc.

# Running
Usage: `./mhdtoy [-o run_name] [-s settings_file] [-S state_file] [-i job_index]`
- `-o` or `--output` allows you to specify the name of the run, which is used as the prefix for the resulting `.out` and `.state` files. Default is `output`.
- `-c` or `--config` allows you to specify the config file to be used for configuring the simulation. Default is `default.config`.
- `-S` or `--state` allows you to specify a .state file from a previous run to be used for initialization. If none specified, a default initialization is used.
- `-i` or `--index` allows you to specify the index of job array index of the current run (for use with Slurm job arrays, see below). Default value is 0.

The resulting simulation will be output as a time series to the file `[run_name].out` and as a state file to `[run_name].state` (which can be used to initialize a future simulation from the end point of that simulation).

# Simulation Settings
All settings are set at run time by the contents of the specified `.config` file. See `default.config` for the required syntax of the `.config` file.