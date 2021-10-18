#include "ucnputils.hpp"

MhdInp gen_inp_grids_ucnp(const PlasmaSettings& pms)
{
    // set up position grids
    Grid x, y, dx, dy;
    int Nx, Ny;
    std::vector<double> dx_vec, dy_vec;
    
    if (strcmp(pms.getopt("grid_opt"),"uniform")){
        Nx = 2*round2int(pms.getvar("x_lim")/pms.getvar("dx"))+1;
        Ny = 2*round2int(pms.getvar("y_lim")/pms.getvar("dy"))+1;
        dx_vec.resize(Nx,pms.getvar("dx"));
        dy_vec.resize(Ny,pms.getvar("dy"));
    }
    else if (strcmp(pms.getopt("grid_opt"),"nonuniform")){
        Nx = 2*floor(pms.getvar("Nx")/2)+1; // ensure number of grids is odd so that there is a 'center' pixel
        Ny = 2*floor(pms.getvar("Ny")/2)+1; // ensure number of grids is odd so that there is a 'center' pixel
        genNonUniformGrids(pms.getvar("x_lim"),pms.getvar("dx"),pms.getvar("Nx"),dx_vec);
        genNonUniformGrids(pms.getvar("y_lim"),pms.getvar("dy"),pms.getvar("Ny"),dy_vec);
    }

    meshgrid(dx_vec,dy_vec,dx,dy);
    std::string center_opt = "center";
    x = PlasmaDomain::convertCellSizesToCellPositions(dx,0,center_opt);
    y = PlasmaDomain::convertCellSizesToCellPositions(dy,1,center_opt);


    MhdInp grids(Nx,Ny);
    double tau_exp { get_tau_exp(pms.getvar("sig_x"),pms.getvar("m_i"),pms.getvar("Te")) };
    grids.set_ion_mass(pms.getvar("m_i"));
    grids.set_adiabatic_index(pms.getvar("gam"));
    grids.set_duration(pms.getvar("tmax/tau_exp")*tau_exp);
    grids.set_var(PlasmaDomain::d_x,dx,center_opt);
    grids.set_var(PlasmaDomain::d_y,dy,center_opt);


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

void genNonUniformGrids(double r_max, double dr_min, int Nr,std::vector<double>& dr)
{
    double p { 10 };
    double f_sum{0};
    double N { ceil(Nr/2) };
    for (int i = 0; i < N; i++){
        f_sum += abs(pow(i,p));
    }
    double m { (r_max - dr_min*N)/f_sum };
    
    dr.resize(Nr,0.);
    for (int i = 0; i < Nr; i++){
        dr[i] = dr_min + abs(m*pow(i-N,p));
    }
}