#ifndef MHD_HPP
#define MHD_HPP

#include "grid.hpp"
#include "plasmadomain.hpp"
#include <filesystem>
namespace fs = std::filesystem;

//Explicit variable specification
void mhdSolve(std::vector<Grid> input_vars, double ion_mass, double adiabatic_index, double time_duration,
              const fs::path &output_pathname, const fs::path &config_filename);

//Specification from complete directory of past run (incl. .state file)
void mhdSolve(const fs::path &prev_run_directory, double time_duration);

#endif