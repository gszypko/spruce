#include "ucnputils.hpp"

MhdInp gen_inp_grids_ucnp(const PlasmaSettings& pms)
{
    Grid x, y;
    int Nx { 2*round2int(pms.getvar("x_lim")/pms.getvar("dx"))+1 };
    int Ny { 2*round2int(pms.getvar("y_lim")/pms.getvar("dy"))+1 };
    std::vector<double> x_vec { linspace<double>(-pms.getvar("x_lim"),pms.getvar("x_lim"),Nx) };
    std::vector<double> y_vec { linspace<double>(-pms.getvar("y_lim"),pms.getvar("y_lim"),Ny) };
    meshgrid(x_vec,y_vec,x,y);

    MhdInp grids(Nx,Ny);
    double tau_exp { get_tau_exp(pms.getvar("sig_x"),pms.getvar("m_i"),pms.getvar("Te")) };
    grids.set_ion_mass(pms.getvar("m_i"));
    grids.set_adiabatic_index(pms.getvar("gam"));
    grids.set_duration(pms.getvar("tmax/tau_exp")*tau_exp);
    grids.set_var(PlasmaDomain::pos_x,x);
    grids.set_var(PlasmaDomain::pos_y,y);
    double rho_max { pms.getvar("n")*pms.getvar("m_i") };
    double rho_min { pms.getvar("n_min")*pms.getvar("m_i") };
    if (strcmp(pms.getopt("n_dist"),"gaussian")){
        grids.set_var(PlasmaDomain::rho,gaussian2D(x,y,rho_max,rho_min,pms.getvar("sig_x"),pms.getvar("sig_y"),0,0));
    }
    else if (strcmp(pms.getopt("n_dist"),"exponential")){
        grids.set_var(PlasmaDomain::rho,exponential2D(x,y,rho_max,rho_min,pms.getvar("sig_x"),pms.getvar("sig_y"),0,0));
    }
    else std::cerr << "<n_dist> must be specified as <gaussian> or <exponential>";

    grids.set_var(PlasmaDomain::temp,Grid(Nx,Ny,pms.getvar("Te")));
    grids.set_var(PlasmaDomain::mom_x,Grid(Nx,Ny,0));
    grids.set_var(PlasmaDomain::mom_y,Grid(Nx,Ny,0));

    std::vector<Grid> B(3.);
    for (int i = 0; i < B.size(); i++){ // for each B-field component
        B[i] = Grid(Nx,Ny); // initialize that component
        for (int j = 0; j < B[i].rows(); j++){
            for (int k = 0; k < B[i].cols(); k++){
                std::vector<double> Bquad = quadrupole_field(x(j,k),y(j,k),0.,pms.getvar("dBdx"));
                B[i](j,k) = Bquad[i];
            }
        }
    } 

    grids.set_var(PlasmaDomain::b_x,B[0]);
    grids.set_var(PlasmaDomain::b_y,B[1]);
    grids.set_var(PlasmaDomain::b_z,B[2]);
    grids.set_var(PlasmaDomain::grav_x,Grid(Nx,Ny,0));
    grids.set_var(PlasmaDomain::grav_y,Grid(Nx,Ny,0));

    return grids;
}