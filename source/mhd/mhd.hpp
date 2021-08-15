#ifndef MHD_HPP
#define MHD_HPP

#include "grid.hpp"
#include "plasmadomain.hpp"
#include <filesystem>

//Explicit variable specification
void mhdSolve(std::vector<Grid> input_vars, const char* output_pathname, const char* config_filename, double time_duration);

//Specification from complete directory of past run (incl. .state file)
void mhdSolve(const char* prev_run_directory, double time_duration);

// //Variable specification from .state file
// void mhdSolve(const char* state_filename, const char* output_pathname, const char* config_filename, double time_duration);

#endif