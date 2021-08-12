#ifndef MHD_HPP
#define MHD_HPP

#include "grid.hpp"

//Explicit variable specification
void mhdSolve(std::vector<Grid> input_vars, const char* output_pathname, const char* config_filename);

//Specification from complete directory of past run (incl. .state file)
void mhdSolve(const char* output_pathname, const char* run_directory);

//Variable specification from .state file
void mhdSolve(const char* state_filename, const char* output_pathname, const char* config_filename);

#endif