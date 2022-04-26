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
    std::string plasma_type = getCommandLineArg(argc, argv, "-t", "--type");

    fs::path prev_run_path(getCommandLineArg(argc, argv, "-p", "--prev"));
    fs::path out_path(getCommandLineArg(argc, argv, "-o", "--output"));
    fs::path config_path(getCommandLineArg(argc, argv, "-c", "--config"));
    fs::path settings_path(getCommandLineArg(argc, argv, "-s", "--settings"));

    int task_array;
    std::string task_array_str = getCommandLineArg(argc, argv, "-i", "--index");
    if(task_array_str.empty()) task_array = 0; //default array index is zero
    else task_array = std::stoi(task_array_str);
    fs::path array_path("array" + std::to_string(task_array));

    double time_duration; //for previous run duration set by command line, for .settings run duration set by .settings
    std::string time_duration_str = getCommandLineArg(argc, argv, "-d", "--duration");
    if(time_duration_str.empty()) time_duration = 1.0; //default duration is 1 second
    else time_duration = std::stod(time_duration_str); //for use with previous run only

    // continue from previous run
    if(!prev_run_path.empty()){
        fs::directory_entry prev_run_dir(prev_run_path);
        assert(prev_run_dir.exists() && prev_run_dir.is_directory() && "Given previous run path must be existing directory");
        mhdSolve(prev_run_path, time_duration);
        return 0;
    }
    // initialize using settings file and pertinent initializer fxn
    else {
        assert(!config_path.empty() && "config file must be specified");
        assert(!settings_path.empty() && "settings file must be specified");
        assert(!out_path.empty() && "output directory must be specified");

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
        time_duration = grids_inp.duration();

        // run MHD
        mhdSolve(grids_inp.grids(), grids_inp.ion_mass(), grids_inp.adiabatic_index(), time_duration, out_path, config_path);
        return 0;
    }
    
}
