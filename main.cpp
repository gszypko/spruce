// Modules: module load ... ... ... -> need to test on cluster eventually, should just need gcc > 8
// Compile: g++ -I source/ -I source/mhd/ -I source/user-interface/ -std=c++17 -fopenmp -g source/ui_main.cpp source/mhd/* source/user-interface/* -o main -lm -lstdc++fs
// Run: run <plasma_type> <out_path> <array_num>
//      plasma_type: "ucnp" or "solar"
//      out_path: relative path to simulation save directory
//      array_num: integer corresponding to row of unique combinations of plasma settings

#include <filesystem>
namespace fs = std::filesystem;
#include <iostream>
#include "mhd.hpp"
#include "plasmadomain.hpp"
#include "utils.hpp"
#include "constants.hpp"

int main(int argc, char *argv[])
{
    // proccess command line arguments
    std::string run_mode = getCommandLineArg(argc, argv, "-m", "--mode");
    std::string time_duration_str = getCommandLineArg(argc, argv, "-d", "--duration");
    double time_duration = time_duration_str.empty() ? -1.0 : std::stod(time_duration_str); //set to -1 if unspecified on cmd line
    std::string cluster_time_str = getCommandLineArg(argc, argv, "-r", "--runtime");
    double cluster_time = cluster_time_str.empty() ? -1.0 : std::stod(cluster_time_str); //set to -1 if unspecified on cmd line
    fs::path out_path(getCommandLineArg(argc, argv, "-o", "--output"));
    assert(!out_path.empty() && "output directory must be specified");

    // handle when running in continue mode
        // time duration must be specified via the command line
        // the previous run directory must exist
    if (run_mode == "continue"){
        fs::path prev_run_path = out_path;
        assert(!time_duration_str.empty() && "In Continue Mode, duration of simulation must be specified on command line");
        fs::directory_entry prev_run_dir(prev_run_path);
        assert(prev_run_dir.exists() && prev_run_dir.is_directory() && "Given output directory of previous run must be existing directory");
        #if VERBOSE 
            std::cout << "Running in Continue Mode for " << time_duration << " s...\n"; 
        #endif
        mhdSolve(prev_run_path, time_duration, cluster_time);
        return 0;
    }
    else if (run_mode == "input"){
        fs::path config_path(getCommandLineArg(argc, argv, "-c", "--config"));
        fs::path grid_path(getCommandLineArg(argc, argv, "-s", "--state"));
        bool seek_config = config_path.empty();
        bool seek_grids = grid_path.empty();
        if(seek_config || seek_grids){
            for(auto const& dir_entry : fs::directory_iterator{out_path}){
                std::string exten = dir_entry.path().extension().string();
                if(seek_config && exten == ".config"){
                    assert(config_path.empty() && "There must be only one .config file in the specified directory");
                    config_path = dir_entry.path();
                    #if VERBOSE
                        std::cout << "Found configuration file " << config_path << std::endl;
                    #endif
                } else if(seek_grids && dir_entry.path().filename() == "init.state"){
                    grid_path = dir_entry.path();
                    #if VERBOSE
                        std::cout << "Found initializing state file " << grid_path << std::endl;
                    #endif
                }
            }
        }
        assert(!config_path.empty() && "Config file not specified and not found in output directory");
        assert(!grid_path.empty() && "Initializing state file not specified and not found in output directory");

        fs::create_directories(out_path);
        if(grid_path.extension().string() == ".state"){ // Specify initial state with .state file
            fs::directory_entry grid_path_dir(grid_path);
            assert(grid_path_dir.exists() && grid_path_dir.is_regular_file() && "Given state file must exist and be a file");
            if(time_duration_str.empty()){
                #if VERBOSE
                std::cout << "Running in Input Mode from the state file " << grid_path.string() << " for duration specified in .config file." << std::endl;
                #endif
            }
            else {
                #if VERBOSE
                std::cout << "Running in Input Mode from the state file " << grid_path.string() << " for " << time_duration << " s...\n";
                #endif
            }
            mhdSolve(grid_path, config_path, out_path, time_duration, !seek_grids, cluster_time);
            return 0;
        } else {
            std::cerr << "Grids must be specified in .state file.\n";
            return 1;
        }
    }
    else {
        std::cerr << "Mode \'" << run_mode << "\' not recognized\n";
        return 1;
    }   
}
