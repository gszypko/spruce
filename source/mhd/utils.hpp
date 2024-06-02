#ifndef UTILS_HPP
#define UTILS_HPP

#include <string>
#include <sstream>
#include <vector>
#include <grid.hpp>

//Erase all occurences of ' ', '\t', and '\n' in str.
//Modifies in-place.
void clearWhitespace(std::string &str);

//Wrapper for using std::getline, then removing all whitespace characters from the result
//Meant to avoid strange behavior when reading directly from a file (e.g. ifstream)
std::istream &getCleanedLine(std::istream &is, std::string &str, char delim = '\n');

//Splits string into vector of substrings, delimited by delim
std::vector<std::string> splitString(const std::string &str, const char delim);

//Retrieve comand line argument preceded by either short_flag or long_flag (both must, at a minimum, start with '-')
std::string getCommandLineArg(int argc, char* argv[], std::string short_flag, std::string long_flag);

double bilinearInterpolate(const std::vector<double> &point, const Grid &quantity, const std::vector<double> &x, const std::vector<double> &y);

#endif
