// Modules: module load ... ... ... -> need to test on cluster eventually, should just need gcc > 8
// Compile: g++ -I source/ -I source/mhd/ -I source/user-interface/ -std=c++17 -fopenmp -g source/ui_main.cpp source/mhd/* source/user-interface/* -o main -lm -lstdc++fs
// Run: run <plasma_type> <out_path> <array_num>
//      plasma_type: "ucnp" or "solar"
//      out_path: relative path to simulation save directory
//      array_num: integer corresponding to row of unique combinations of plasma settings

#include "MhdInp.hpp"
#include "ui_utility.hpp"
#include <filesystem>
namespace fs = std::filesystem;
#include "PlasmaSettings.hpp"
#include "mhd.hpp"
#include "utils.hpp"
#include "plasmadomain.hpp"
#include "solarutils.hpp"
#include "ucnputils.hpp"

int main(int argc, char *argv[])
{
    // unfold executable input parameters
    // TODO: encapsulate command line arg parsing into class (make checking for correct specification easier)
    std::string time_duration_str = getCommandLineArg(argc, argv, "-d", "--duration");
    double time_duration = time_duration_str.empty() ? -1.0 : std::stod(time_duration_str); //set to -1 if unspecified on cmd line

    // Continue Mode
    fs::path prev_run_path(getCommandLineArg(argc, argv, "-p", "--prev"));
    if(!prev_run_path.empty()){
        assert(!time_duration_str.empty() && "In Continue Mode, duration of simulation must be specified on command line");
        fs::directory_entry prev_run_dir(prev_run_path);
        assert(prev_run_dir.exists() && prev_run_dir.is_directory() && "Given previous run path must be existing directory");

        std::cout << "Running in Continue Mode for " << time_duration << " s...\n";
        mhdSolve(prev_run_path, time_duration);
        return 0;
    }

    fs::path out_path(getCommandLineArg(argc, argv, "-o", "--output"));
    fs::path config_path(getCommandLineArg(argc, argv, "-c", "--config"));
    assert(!config_path.empty() && "config file must be specified");
    assert(!out_path.empty() && "output directory must be specified");

    // Custom Input Mode
    fs::path grid_path(getCommandLineArg(argc, argv, "-g", "--grids"));
    if(!grid_path.empty()){
        fs::create_directories(out_path);
        if(grid_path.extension().string() == ".state"){ // Specify initial state with .state file
            fs::directory_entry grid_path_dir(grid_path);
            assert(grid_path_dir.exists() && grid_path_dir.is_regular_file() && "Given state file must exist and be a file");
            if(time_duration_str.empty()){
                std::cout << "Running in Custom Input Mode from the state file " << grid_path.string() << " for duration specified in file...\n";
            }
            else {
                std::cout << "Running in Custom Input Mode from the state file " << grid_path.string() << " for " << time_duration << " s...\n";
            }
            mhdSolve(grid_path, config_path, out_path, time_duration);
            return 0;
        } else {
            std::cerr << "Grids must be specified in .state file.\n";
            return 1;
        }
    }

    // Problem Generator Mode
    std::string plasma_type = getCommandLineArg(argc, argv, "-t", "--type");
    if(!plasma_type.empty()){
        fs::path settings_path(getCommandLineArg(argc, argv, "-s", "--settings"));
        assert(!settings_path.empty() && "settings file must be specified");

        int task_array;
        std::string task_array_str = getCommandLineArg(argc, argv, "-i", "--index");
        if(task_array_str.empty()) task_array = 0; //default array index is zero
        else task_array = std::stoi(task_array_str);
        fs::path array_path("array" + std::to_string(task_array));

        // create subdirectory for array index, if one was explicitly given
        if(!task_array_str.empty()) out_path/=array_path;
        fs::create_directories(out_path);

        // load .settings file
        PlasmaSettings pms(settings_path,task_array);
        pms.write_array_params(out_path);
        
        // generate input grids for MHD code
        MhdInp grids_inp;

        if (plasma_type.compare("ucnp") == 0) grids_inp = gen_inp_grids_ucnp(pms);
        else if (plasma_type.compare("solar") == 0) grids_inp = SolarUtils::SolarMHDInput(pms);
        else std::cerr << "Plasma type not specified ('ucnp' or 'solar')";
        if(time_duration_str.empty()) time_duration = grids_inp.duration();

        // run 
        std::cout << "Running in Problem Generator Mode using the " << plasma_type << " problem generator for " << time_duration << " s...\n";
        mhdSolve(grids_inp.grids(), grids_inp.ion_mass(), grids_inp.adiabatic_index(), time_duration, out_path, config_path);
        return 0;
    }
    
}
