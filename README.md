# mhdtoy
A two-dimensional quasi-MHD simulation for static magnetic field

# Compiling
The Makefile is very dumb and doesn't keep track of any dependencies. It only has phony targets clang, intel, and gnu corresponding to full compilation commands using clang++, icc, and gcc compilers respectively. For instance, execute `make intel` to compile using icc.

# Running
Executing the resulting executable (`mhdtoy` or `mhdtoy.exe`) with zero command line arguments will initialize the system with an initial hydrostatic solar atmosphere. The resulting simulation will be output as a time series to the file `output.out` and as a state file to `output.state` (which can be used to initialize a future simulation from the end point of that simulation; see below).

Executing with a single command line argument will do the same as above, but with the given argument being used for the resulting file names. For instance, `./mhdtoy testing` will give files `testing.out` and `testing.state`.

Executing with two command line arguments allows you to specify a `.state` file to initialize the simulation. The first argument specifies the resulting `.out` and `.state` file name (as above) and the second argument is the filename of the `.state` file to be used to initialize the system.

# Simulation Settings
Currently, all of the settings for the simulation are specified in the `constants.hpp` file and are therefore specified at compile time. (Currently working toward specifying settings at execution time, to avoid the need for multiple executables.)
