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

//Erase all occurences of ' ', '\t', and '\n' in str.
//Modifies in-place.
void clearWhitespace(std::string &str);

//Splits string into vector of substrings, delimited by delim
std::vector<std::string> splitString(const std::string &str, const char delim);

//Retrieve comand line argument preceded by either short_flag or long_flag (both must, at a minimum, start with '-')
std::string getCommandLineArg(int argc, char* argv[], std::string short_flag, std::string long_flag);

#endif
