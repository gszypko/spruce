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

MhdInp gen_inp_grids_ucnp(const PlasmaSettings& pms)
{
    Grid x, y;
    int Nx { round2int(2*pms.getvar("x_lim")/pms.getvar("dx")) };
    int Ny { round2int(2*pms.getvar("y_lim")/pms.getvar("dy")) };
    std::vector<double> x_vec { linspace<double>(-pms.getvar("x_lim"),pms.getvar("x_lim"),Nx) };
    std::vector<double> y_vec { linspace<double>(-pms.getvar("y_lim"),pms.getvar("y_lim"),Ny) };
    meshgrid(x_vec,y_vec,x,y);

    MhdInp grids(Nx,Ny);
    grids.set_var(PlasmaDomain::pos_x,x);
    grids.set_var(PlasmaDomain::pos_y,y);
    if (strcmp(pms.getopt("n_dist"),"gaussian")){
        grids.set_var(PlasmaDomain::rho,gaussian2D(x,y,pms.getvar("n"),pms.getvar("sig_x"),pms.getvar("sig_y"),0,0));
    }
    else if (strcmp(pms.getopt("n_dist"),"exponential")){
        grids.set_var(PlasmaDomain::rho,exponential2D(x,y,pms.getvar("n"),pms.getvar("sig_x"),pms.getvar("sig_y"),0,0));
    }
    else std::cerr << "<n_dist> must be specified as <gaussian> or <exponential>";

    grids.set_var(PlasmaDomain::temp,Grid(Nx,Ny,pms.getvar("Te")));
    grids.set_var(PlasmaDomain::mom_x,Grid(Nx,Ny,0));
    grids.set_var(PlasmaDomain::mom_y,Grid(Nx,Ny,0));

    std::vector<Grid> B(3.);
    for (int i = 0; i < B.size(); i++){
        B[i] = Grid(Nx,Ny);
        for (int j = 0; j < B[i].rows(); j++){
            for (int k = 0; k < B[i].cols(); k++){
                std::vector<double> Bquad = quadrupole_field(x(j,k),y(j,k),0.,pms.getvar("dBdx"));
                B[i](j,k) = Bquad[i];
            }
        }
    } 

    grids.set_var(PlasmaDomain::b_x,B[0]);
    grids.set_var(PlasmaDomain::b_y,B[1]);
    grids.set_var(PlasmaDomain::b_z,Grid(Nx,Ny,0));

    return grids;
}

int main(int argc, char *argv[])
{
    // unfold executable input parameters
    std::filesystem::path prev_run_path = getCommandLineArg(argc, argv, "-p", "--prev");

    std::string plasma_type = getCommandLineArg(argc, argv, "-t", "--type");

    std::filesystem::path out_path = getCommandLineArg(argc, argv, "-o", "--output");
 
    std::string config_path = getCommandLineArg(argc, argv, "-c", "--config");

    int task_array;
    std::string task_array_str = getCommandLineArg(argc, argv, "-i", "--index");
    if(task_array_str.empty()) task_array = 0; //default array index is zero
    else task_array = std::stoi(task_array_str);

    double time_duration;
    std::string time_duration_str = getCommandLineArg(argc, argv, "-d", "--duration");
    if(time_duration_str.empty()) time_duration = 1.0; //default duration is 1 second
    else time_duration = std::stod(time_duration_str);


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
        if(!task_array_str.empty()) out_path /= ("array"+std::to_string(task_array));
        std::filesystem::create_directories(out_path);

        // choose which plasma settings file to use
        std::filesystem::path settings_path{}; // relative path to .settings file
        if (plasma_type.compare("ucnp") == 0) settings_path = "ucnp.settings";
        else if (plasma_type.compare("solar") == 0) settings_path = "solar.settings";
        else std::cerr << "Plasma type not specified ('ucnp' or 'solar')";

        // load .settings file
        PlasmaSettings pms(settings_path,task_array);
        pms.write_array_params(out_path);
        
        // generate input grids for MHD code
        MhdInp grids_inp;
        if (plasma_type.compare("ucnp") == 0) grids_inp = gen_inp_grids_ucnp(pms);
        else if (plasma_type.compare("solar") == 0) grids_inp = SolarUtils::SolarMHDInput(pms);
        else std::cerr << "Plasma type not specified ('ucnp' or 'solar')";

        // run MHD
        mhdSolve(grids_inp.grids(), grids_inp.ion_mass(), grids_inp.adiabatic_index(), time_duration, out_path.c_str(), config_path.c_str());
        return 0;
    }
    
}
