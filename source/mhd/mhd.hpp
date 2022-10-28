#ifndef MHD_HPP
#define MHD_HPP

#include "grid.hpp"
#include <filesystem>
namespace fs = std::filesystem;
#include <vector>

//Problem Generator Mode
void mhdSolve(std::vector<Grid> input_vars, double ion_mass, double adiabatic_index, double time_duration,
              const fs::path &output_pathname, const fs::path &config_filename);

//Continue Mode
void mhdSolve(const fs::path &prev_run_directory, double time_duration, double cluster_time);

//Custom Input Mode (from state file)
void mhdSolve(const fs::path &state_filename, const fs::path &config_filename, const fs::path &output_pathname, double time_duration, bool overwrite_init, double cluster_time);

#endif