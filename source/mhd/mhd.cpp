#include <filesystem>
#include "grid.hpp"
#include "plasmadomain.hpp"

void mhd(const Grid& rho, const Grid& temp, const Grid& mom_x, const Grid& mom_y,
         const Grid& b_x, const Grid& b_y, const Grid& b_z, const Grid& pos_x, const Grid& pos_y,
         const Grid& grav_x, const Grid& grav_y, const char* output_pathname, const char* config_filename)
{
  PlasmaDomain simulation(output_pathname,config_filename);
  simulation.initialize(rho, temp, mom_x, mom_y, b_x, b_y, b_z, pos_x, pos_y, grav_x, grav_y);
  simulation.run();
}

void mhd(const char* state_filename, const char* output_pathname, const char* config_filename)
{
  PlasmaDomain simulation(output_pathname,config_filename);
  simulation.readStateFile(state_filename);
  simulation.run();
}

void mhd(const char* output_pathname, const char* run_directory)
{
  std::string state_filename, config_filename;
  std::filesystem::path run_dir_path(run_directory);

  for(auto const& dir_entry : std::filesystem::directory_iterator{run_dir_path}){
    std::string exten = dir_entry.path().extension();
    if(exten == ".config"){
      assert(config_filename.empty() && "There must be only one .config file in the specified directory");
      config_filename = dir_entry.path();
    } else if(exten == ".state"){
      assert(state_filename.empty() && "There must be only one .state file in the specified directory");
      state_filename = dir_entry.path();
    }
  }

  assert(!config_filename.empty() && !state_filename.empty() && "There must be a .state and a .config file in the specified directory");

  mhd(state_filename.c_str(), output_pathname, config_filename.c_str());
}