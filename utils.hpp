#ifndef UTILS_HPP
#define UTILS_HPP

#include "grid.hpp"

//Generates gaussian initial condition for a variable, centered at middle of grid
Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y);

//Generates potential bipolar field for component corresponding to index "index"
//Centered s.t. origin lies at bottom middle of domain
//Pressure scale height h, field poles at +/- l, field strength at poles b0
Grid BipolarField(int xdim, int ydim, double b0, double h, double dx, double dy, int index);

//Generates grid with exponential falloff in the y-direction, with the quantity
//"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
Grid HydrostaticFalloff(double base_value, double scale_height, int xdim, int ydim, double dy);

//Erase all occurences of ' ', '\t', and '\n' in str.
//Modifies in-place.
void clearWhitespace(std::string &str);

//Splits string into vector of substrings, delimited by delim
std::vector<std::string> splitString(const std::string &str, const char delim);

#endif
