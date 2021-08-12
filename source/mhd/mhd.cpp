#include <filesystem>
#include "grid.hpp"
#include "plasmadomain.hpp"

void mhdSolve(std::vector<Grid> input_vars, const char* output_pathname, const char* config_filename)
{
  PlasmaDomain simulation(output_pathname,config_filename);
  simulation.initialize(input_vars);
  simulation.run();
}

void mhdSolve(const char* state_filename, const char* output_pathname, const char* config_filename)
{
  PlasmaDomain simulation(output_pathname,config_filename);
  simulation.readStateFile(state_filename);
  simulation.run();
}

void mhdSolve(const char* output_pathname, const char* run_directory)
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

  mhdSolve(state_filename.c_str(), output_pathname, config_filename.c_str());
}