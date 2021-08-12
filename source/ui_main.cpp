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

MhdInp gen_inp_grids_ucnp(const PlasmaSettings& pms)
{
    Grid x, y;
    int Nx { round2int(2*pms.getvar("x_lim")/pms.getvar("dx")) };
    int Ny { round2int(2*pms.getvar("y_lim")/pms.getvar("dy")) };
    std::vector<double> x_vec { linspace<double>(-pms.getvar("x_lim"),pms.getvar("x_lim"),Nx) };
    std::vector<double> y_vec { linspace<double>(-pms.getvar("y_lim"),pms.getvar("y_lim"),Ny) };
    meshgrid(x_vec,y_vec,x,y);

    MhdInp grids(Nx,Ny);
    grids.set_var(MhdInp::x,x);
    grids.set_var(MhdInp::y,y);
    if (strcmp(pms.getopt("n_dist"),"gaussian")){
        grids.set_var(MhdInp::rho,gaussian2D(x,y,pms.getvar("n"),pms.getvar("sig_x"),pms.getvar("sig_y"),0,0));
    }
    else if (strcmp(pms.getopt("n_dist"),"exponential")){
        grids.set_var(MhdInp::rho,exponential2D(x,y,pms.getvar("n"),pms.getvar("sig_x"),pms.getvar("sig_y"),0,0));
    }
    else std::cerr << "<n_dist> must be specified as <gaussian> or <exponential>";

    grids.set_var(MhdInp::temp,Grid(Nx,Ny,pms.getvar("Te")));
    grids.set_var(MhdInp::mom_x,Grid(Nx,Ny,0));
    grids.set_var(MhdInp::mom_y,Grid(Nx,Ny,0));

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

    grids.set_var(MhdInp::b_x,B[0]);
    grids.set_var(MhdInp::b_y,B[1]);
    grids.set_var(MhdInp::b_z,Grid(Nx,Ny,0));

    return grids;
}

int main(int argc, char *argv[])
{
    // unfold executable input parameters
    std::string plasma_type = argv[1]; // .settings option - "ucnp" or "solar"
    std::filesystem::path out_path = argv[2]; // relative save directory path
    int task_array{atoi(argv[3])};  // array number from slurm file input
    std::string config_path = argv[4]; // path of .config file to use

    //refactor mhd function to take MhdInp, out_path, config_path
    //command line flag parsing? so that you can either specify all of the above or specify a previous run directory

    // create subdirectory for array index
    out_path /= ("array"+std::to_string(task_array));
    std::filesystem::create_directories(out_path);

    // choose which plasma settings file to use
    std::filesystem::path settings_path{}; // relative path to .settings file
    if (plasma_type.compare("ucnp") == 0) settings_path = "ucnp.settings";
    else if (plasma_type.compare("solar") == 0) settings_path = "solar.settings";
    else std::cerr << "First argument to executable must be <ucnp> or <solar>";

    // load .settings file
    PlasmaSettings pms(settings_path,task_array);
    pms.write_array_params(out_path);
    
    // generate input grids for MHD code - if statement to select plasma type
    MhdInp grids_inp = gen_inp_grids_ucnp(pms);
    // grids_inp.all_initialized();

    // run MHD
    mhdSolve(grids_inp.grids(),out_path.c_str(),config_path.c_str());

    return 0;
}
