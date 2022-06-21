#include "ucnputils.hpp"
MhdInp gen_inp_grids_ucnp(const Settings& pms)
{
    // set up position grids
    Grid x, y, dx, dy;
    int Nx, Ny;
    std::vector<double> dx_vec, dy_vec;
    
    Nx = 2*floor(pms.getval("Nx")/2)+1; // ensure number of grids is odd so that there is a 'center' pixel
    Ny = 2*floor(pms.getval("Ny")/2)+1; // ensure number of grids is odd so that there is a 'center' pixel
    genNonUniformGrids(pms.getval("x_lim"),Nx,dx_vec,pms.getopt("grid_opt"));
    genNonUniformGrids(pms.getval("y_lim"),Ny,dy_vec,pms.getopt("grid_opt"));

    meshgrid(dx_vec,dy_vec,dx,dy);
    std::string center_opt = "center";
    x = PlasmaDomain::convertCellSizesToCellPositions(dx,0,center_opt);
    y = PlasmaDomain::convertCellSizesToCellPositions(dy,1,center_opt);

    PlasmaDomain pd {};
    MhdInp grids(Nx,Ny,pd,"ideal_mhd");
    double sig = pow(pms.getval("sig_x")*pow(pms.getval("sig_y"),2.),1./3.); // geometric mean of plasma size
    double tau { tau_exp(sig,pms.getval("m_i"),2*pms.getval("Te")) };
    grids.set_ion_mass(pms.getval("m_i"));
    grids.set_adiabatic_index(pms.getval("gam"));
    grids.set_duration(pms.getval("t_max/tau_exp")*tau);
    grids.set_time_output_interval(pms.getval("t_iter/tau_exp")*tau);
    grids.set_var("d_x",dx);
    grids.set_var("d_y",dy);
    grids.set_var("pos_x",x);
    grids.set_var("pos_y",y);

    double rho_max { pms.getval("n")*pms.getval("m_i") };
    double rho_min { pms.getval("n_min")*pms.getval("m_i") };
    if (strcmp(pms.getopt("n_dist"),"gaussian")){
        grids.set_var("rho",gaussian2D(x,y,rho_max,rho_min,pms.getval("sig_x"),pms.getval("sig_y"),0,0));
    }
    else if (strcmp(pms.getopt("n_dist"),"exponential")){
        grids.set_var("rho",exponential2D(x,y,rho_max,rho_min,pms.getval("sig_x"),pms.getval("sig_y"),0,0));
    }
    else std::cerr << "<n_dist> must be specified as <gaussian> or <exponential>";

    grids.set_var("temp",Grid(Nx,Ny,pms.getval("Te")));
    grids.set_var("mom_x",Grid(Nx,Ny,0));
    grids.set_var("mom_y",Grid(Nx,Ny,0));

    // initialize quadrupole magnetic fields
    AntiHelmholtz quad(30.,120.,pms.getval("dBdx"),CurrentLoop::ax_x);

    std::vector<Grid> B(2.);
    for (int i = 0; i < B.size(); i++){ // for each B-field component (x and y)
        B[i] = Grid::Zero(Nx,Ny); // initialize that component
        for (int j = 0; j < B[i].rows(); j++){
            for (int k = 0; k < B[i].cols(); k++){
                std::vector<long double> Bquad = quad.get_field({x(j,k),y(j,k),0.});
                B[i](j,k) = Bquad[i];
            }
        }
    } 

    grids.set_var("be_x",B[0]);
    grids.set_var("be_y",B[1]);
    grids.set_var("bi_x",Grid(Nx,Ny,0));
    grids.set_var("bi_y",Grid(Nx,Ny,0));
    grids.set_var("grav_x",Grid(Nx,Ny,0));
    grids.set_var("grav_y",Grid(Nx,Ny,0));

    return grids;
}

void genNonUniformGrids(double r_max, int Nr,std::vector<double>& dr,std::string opt)
{
    double A = 1.04;
    double B;
    if (strcmp(opt,"uniform")) B = 1;
    else if (strcmp(opt,"nonuniform")) B = 0.995;
    double N = ceil(Nr/2.);
    std::vector<double> dx(N);
    std::vector<double> x(N);
    dx[0] = 1;
    x[0] = 0;
    for (int i = 1; i < N; i++){
        dx[i] = 1 + A*(dx[i-1]-B);
        x[i] = x[i-1] + dx[i]/2. + dx[i-1]/2.;
    }
    double x_max = x.back();
    
    dr.clear();
    dr.reserve(Nr);
    for (int i = N-1; i > 0; i--){
        dr.push_back(dx[i]/x_max*r_max);
    }
    for (int i = 0; i < N; i++){
        dr.push_back(dx[i]/x_max*r_max);
    }
}
