#include "mhd.hpp"
#include "plasmadomain.hpp"
#include <cassert>
#include <string>

void mhdSolve(const fs::path &prev_run_directory, double time_duration)
{
  fs::path state_filename, config_filename;

  for(auto const& dir_entry : fs::directory_iterator{prev_run_directory}){
    std::string exten = dir_entry.path().extension().string();
    if(exten == ".config"){
      assert(config_filename.empty() && "There must be only one .config file in the specified directory");
      config_filename = dir_entry.path();
    } else if(dir_entry.path().filename() == "end.state"){
      assert(state_filename.empty() && "The state file to continue from must be named end.state");
      state_filename = dir_entry.path();
    }
  }

  assert(!config_filename.empty() && !state_filename.empty() && "There must be a .state and a .config file in the specified directory");

  PlasmaDomain simulation(prev_run_directory,config_filename,state_filename,true,false); // run in continue mode
  simulation.run(time_duration);
}

void mhdSolve(const fs::path &state_filename, const fs::path &config_filename, const fs::path &output_pathname, double time_duration, bool overwrite_init)
{
  PlasmaDomain simulation(output_pathname,config_filename,state_filename,false,overwrite_init); // run from state file input
  simulation.run(time_duration);
}
