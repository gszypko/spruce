// Modules: module load ... ... ... -> need to test on cluster eventually, should just need gcc > 8
// Compile: g++-8 -I source/ -I source/mhd/ -I source/user-interface/ -std=c++17 -fopenmp -g source/ui_main.cpp source/mhd/* source/user-interface/* -o run -lm -lstdc++fs
// Run: run <plasma_type> <out_path> <array_num>
//      plasma_type: "ucnp" or "solar"
//      out_path: relative path to simulation save directory
//      array_num: integer corresponding to row of unique combinations of plasma settings

#include "MhdInp.hpp"
#include "ui_utility.hpp"
#include <filesystem>
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
    std::filesystem::path prev_run_path = getCommandLineArg(argc, argv, "-p", "--prev");

    std::string plasma_type = getCommandLineArg(argc, argv, "-t", "--type");

    std::filesystem::path out_path = getCommandLineArg(argc, argv, "-o", "--output");
 
    std::string config_path = getCommandLineArg(argc, argv, "-c", "--config");
    if(config_path.empty()) config_path = "default.config";

    std::filesystem::path settings_path = getCommandLineArg(argc, argv, "-s", "--settings");

    int task_array;
    std::string task_array_str = getCommandLineArg(argc, argv, "-i", "--index");
    if(task_array_str.empty()) task_array = 0; //default array index is zero
    else task_array = std::stoi(task_array_str);
    std::string array_path = "array" + std::to_string(task_array);

    double time_duration; //for previous run duration set by command line, for .settings run duration set by .settings
    std::string time_duration_str = getCommandLineArg(argc, argv, "-d", "--duration");
    if(time_duration_str.empty()) time_duration = 1.0; //default duration is 1 second
    else time_duration = std::stod(time_duration_str); //for use with previous run only

    // continue from previous run
    if(!prev_run_path.empty()){
        std::filesystem::directory_entry prev_run_dir(prev_run_path);
        assert(prev_run_dir.exists() && prev_run_dir.is_directory() && "Given previous run path must be existing directory");
        mhdSolve(prev_run_path.c_str(), time_duration);
        return 0;
    }
    // initialize using settings file and pertinent initializer fxn
    else {
        // create subdirectory for array index, if one was explicitly given
        if(!task_array_str.empty()) out_path/=array_path;
        std::filesystem::create_directories(out_path);

        // choose which plasma settings file to use as a default, if none specified
        if (plasma_type.compare("ucnp") == 0){
            if(settings_path.string().empty()){
                settings_path = "ucnp.settings";
            }
        }
        else if (plasma_type.compare("solar") == 0){
            if(settings_path.string().empty()){
                settings_path = "solar.settings";
            }
        }  
        else std::cerr << "Plasma type not specified ('ucnp' or 'solar')\n";

        // load .settings file
        PlasmaSettings pms(settings_path,task_array);
        pms.write_array_params(out_path);
        
        // generate input grids for MHD code
        MhdInp grids_inp;
        if (plasma_type.compare("ucnp") == 0) grids_inp = gen_inp_grids_ucnp(pms);
        else if (plasma_type.compare("solar") == 0) grids_inp = SolarUtils::SolarMHDInput(pms);
        else std::cerr << "Plasma type not specified ('ucnp' or 'solar')";
        time_duration = grids_inp.duration();

        // run MHD
        mhdSolve(grids_inp.grids(), grids_inp.ion_mass(), grids_inp.adiabatic_index(), time_duration, out_path.c_str(), config_path.c_str());
        return 0;
    }
    
}
