#include "mhd.hpp"

void mhdSolve(std::vector<Grid> input_vars, double ion_mass, double adiabatic_index, double time_duration, const char* output_pathname, const char* config_filename)
{
  PlasmaDomain simulation(output_pathname,config_filename);
  simulation.initialize(input_vars,ion_mass,adiabatic_index);
  simulation.run(time_duration);
}

void mhdSolve(const char* prev_run_directory, double time_duration)
{
  std::string state_filename, config_filename;
  std::filesystem::path run_dir_path(prev_run_directory);

  for(auto const& dir_entry : std::filesystem::directory_iterator{run_dir_path}){
    std::string exten = dir_entry.path().extension();
    if(exten == ".config"){
      assert(config_filename.empty() && "There must be only one .config file in the specified directory");
      config_filename = dir_entry.path();
    } else if(dir_entry.path().filename() == "mhd.state"){
      assert(state_filename.empty() && "There must be only one .state file in the specified directory");
      state_filename = dir_entry.path();
    }
  }

  assert(!config_filename.empty() && !state_filename.empty() && "There must be a .state and a .config file in the specified directory");

  PlasmaDomain simulation(prev_run_directory,config_filename.c_str(), true); // run in continue mode
  simulation.readStateFile(state_filename.c_str());
  simulation.run(time_duration);
}