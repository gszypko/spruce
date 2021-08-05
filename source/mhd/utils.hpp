#ifndef UTILS_HPP
#define UTILS_HPP

#include "grid.hpp"

//Generates gaussian initial condition for a variable, centered at middle of grid
Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y);

// //Generates potential bipolar field for component corresponding to index "index"
// //Centered s.t. origin lies at bottom middle of domain
// //Pressure scale height h, field poles at +/- l, field strength at poles b0
// Grid BipolarField(const Grid& m_pos_x, const Grid& m_pos_y, double b0, double h, int index);

// //Generates grid with exponential falloff in the y-direction, with the quantity
// //"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
// Grid HydrostaticFalloff(double base_value, double scale_height, const Grid& m_pos_y);

//Erase all occurences of ' ', '\t', and '\n' in str.
//Modifies in-place.
void clearWhitespace(std::string &str);

//Splits string into vector of substrings, delimited by delim
std::vector<std::string> splitString(const std::string &str, const char delim);

//Returns vector<string> containing the name of the run, the name of the config file,
//and the name of the state file (in that order).
//Run name must be preceded by -o or --out, config file must be preceded by -c or --config,
//state file name must be preceded by -s or --state.
//Default run name is "output" and default config file is "default.config".
//If no state file provided, that entry is left as an empty string.
std::vector<std::string> parseCommandLineArgs(int argc, char* argv[]);

#endif
