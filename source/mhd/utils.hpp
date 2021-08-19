#ifndef UTILS_HPP
#define UTILS_HPP

#include <omp.h>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <algorithm>
#include <sstream>
#include "grid.hpp"
#include "constants.hpp"

//Generates gaussian initial condition for a variable, centered at middle of grid
Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y);

//Erase all occurences of ' ', '\t', and '\n' in str.
//Modifies in-place.
void clearWhitespace(std::string &str);

//Splits string into vector of substrings, delimited by delim
std::vector<std::string> splitString(const std::string &str, const char delim);

// //Returns vector<string> containing the name of the run, the name of the config file,
// //and the name of the state file (in that order).
// //Run name must be preceded by -o or --out, config file must be preceded by -c or --config,
// //state file name must be preceded by -s or --state.
// //Default run name is "output" and default config file is "default.config".
// //If no state file provided, that entry is left as an empty string.
// std::vector<std::string> parseCommandLineArgs(int argc, char* argv[]);

//Retrieve comand line argument preceded by either short_flag or long_flag (both must, at a minimum, start with '-')
std::string getCommandLineArg(int argc, char* argv[], std::string short_flag, std::string long_flag);

#endif
