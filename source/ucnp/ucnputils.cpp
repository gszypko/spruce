#include "ucnputils.hpp"
MhdInp gen_inp_grids_ucnp(const std::unique_ptr<Settings>& pms)
{
    // set up position grids
    Grid x, y, dx, dy;
    int Nx, Ny;
    std::vector<double> dx_vec, dy_vec;
    
    Nx = 2*floor(pms->getval("Nx")/2)+1; // ensure number of grids is odd so that there is a 'center' pixel
    Ny = 2*floor(pms->getval("Ny")/2)+1; // ensure number of grids is odd so that there is a 'center' pixel
    genNonUniformGrids(pms->getval("x_lim"),Nx,dx_vec,pms->getopt("grid_opt"));
    genNonUniformGrids(pms->getval("y_lim"),Ny,dy_vec,pms->getopt("grid_opt"));

    meshgrid(dx_vec,dy_vec,dx,dy);
    std::string center_opt = "center";
    x = PlasmaDomain::convertCellSizesToCellPositions(dx,0,center_opt);
    y = PlasmaDomain::convertCellSizesToCellPositions(dy,1,center_opt);

    PlasmaDomain pd {};
    MhdInp grids(Nx,Ny,pd,"ideal_mhd");
    grids.set_ion_mass(pms->getval("m_i"));
    grids.set_adiabatic_index(pms->getval("gam"));
    grids.set_duration(pms->getval("duration"));
    grids.set_time_output_interval(pms->getval("time_output_interval"));
    grids.set_var("d_x",dx);
    grids.set_var("d_y",dy);
    grids.set_var("pos_x",x);
    grids.set_var("pos_y",y);

    double rho_max { pms->getval("n")*pms->getval("m_i") };
    double rho_min { pms->getval("n_min")*pms->getval("m_i") };
    if (pms->getopt("n_dist")=="gaussian"){
        grids.set_var("rho",gaussian2D(x,y,rho_max,rho_min,pms->getval("sig_x"),pms->getval("sig_y"),0,0));
    }
    else if (pms->getopt("n_dist")=="exponential"){
        grids.set_var("rho",exponential2D(x,y,rho_max,rho_min,pms->getval("sig_x"),pms->getval("sig_y"),0,0));
    }
    else std::cerr << "<n_dist> must be specified as <gaussian> or <exponential>";

    grids.set_var("temp",Grid(Nx,Ny,pms->getval("Te")));
    grids.set_var("mom_x",Grid(Nx,Ny,0));
    grids.set_var("mom_y",Grid(Nx,Ny,0));

    // initialize quadrupole magnetic fields
    AntiHelmholtz quad(30.,120.,pms->getval("dBdx"),CurrentLoop::ax_x);

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
    if (opt=="uniform") B = 1;
    else if (opt=="nonuniform") B = 0.995;
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

// return <x> and <y> mesh grids 
void meshgrid(const std::vector<double>& v_x,const std::vector<double>& v_y,Grid& grid_x,Grid& grid_y)
{
    grid_x = Grid(v_x.size(),v_y.size(),0.);
    grid_y = Grid(v_x.size(),v_y.size(),0.);

    for (int i = 0; i < v_x.size(); i++){
        for (int j = 0; j < v_y.size(); j++){
            grid_x(i,j) = v_x[i];
            grid_y(i,j) = v_y[j];
        }
    }
}

// define 2D gaussian distribution - spatial units must match input grid, output grid units match <amp>
Grid gaussian2D(Grid x,Grid y,double amp,double min,double sig_x,double sig_y,double x_cen,double y_cen)
// x: input meshgrid of x-positions
// y: input meshgrid of y-positions
// amp: amplitude of 2D Gaussian, can be +/-
// sig_x: RMS width of Gaussian distribution on x-axis
// sig_y: RMS width of Gaussian distribution on y-axis
// x_cen: center of Gaussian distribution on x-axis
// y_cen: center of Gaussian distribution on y-axis
// All position units must be the same, otherwise they don't matter. Output grid has same units as <amp>.
{
    bool size_check { true };
    if (x.rows() != y.rows()) size_check = false;
    if (x.cols() != y.cols()) size_check = false;
    if (x.size() != y.size()) size_check = false;
    if (!size_check) std::cerr << "Input grids <x> and <y> must have the same size." << std::endl;
    Grid g2D { Grid(x.rows(),x.cols()) };
    for (int i = 0; i < g2D.rows(); i++){
        for (int j = 0; j < g2D.cols(); j++){
            g2D(i,j) = min+amp*exp(-0.5*pow((x(i,j)-x_cen)/abs(sig_x),2))*exp(-0.5*pow((y(i,j)-y_cen)/abs(sig_y),2));
        }
    }
    return g2D;
}

// define 2D exponential distribution - spatial units must match input grid, output grid units match <amp>
Grid exponential2D(Grid x,Grid y,double amp,double min,double sig_x,double sig_y,double x_cen,double y_cen)
// x: input meshgrid of x-positions
// y: input meshgrid of y-positions
// amp: amplitude of 2D Gaussian, can be +/-
// sig_x: RMS width of Gaussian distribution on x-axis
// sig_y: RMS width of Gaussian distribution on y-axis
// x_cen: center of Gaussian distribution on x-axis
// y_cen: center of Gaussian distribution on y-axis
// All position units must be the same, otherwise they don't matter. Output grid has same units as <amp>.
{
    bool size_check { true };
    if (x.rows() != y.rows()) size_check = false;
    if (x.cols() != y.cols()) size_check = false;
    if (x.size() != y.size()) size_check = false;
    if (!size_check) std::cerr << "Input grids <x> and <y> must have the same size." << std::endl;
    Grid exp2D { Grid(x.rows(),x.cols()) };
    for (int i = 0; i < exp2D.rows(); i++){
        for (int j = 0; j < exp2D.cols(); j++){
            exp2D(i,j) = min+amp*exp(-abs(x(i,j)-x_cen)/abs(sig_x))*exp(-abs(y(i,j)-y_cen)/abs(sig_y));
        }
    }
    return exp2D;
}

namespace phys
{
// ratio between average Coulomb interaction energy and thermal energy (inputs in cgs)
double coulomb_coupling(double n,double T)
// n: plasma density (cm^-3)
// T: temperature (K) of a particular species (i.e., electron or ion)
{ return nearest_coulomb_pot(n)/(K_B*T); }

// average interparticle spacing of gas with density n (SI units)
double wigner_seitz_radius(double n) 
// n: plasma density (cm^-3)
{ return pow(3./(4.*PI*n),1./3.); } 

// plasma oscillation frequency (rad/s, inputs cgs)
double plasma_freq(double n,double m) 
// n: plasma density (cm^-3)
// m: particle mass (g)
{ return sqrt(4*PI*n*pow(E,2.)/m); }

// einsten frequency (rad/s, inputs cgs)
double einstein_freq(double n,double m)
// n: plasma density (cm^-3)
// m: particle mass (g)
{ return plasma_freq(n,m)/sqrt(3.); }

// average Coulomb potential between nearest neighbors, cgs units
double nearest_coulomb_pot(double n) 
// n: plasma density (cm^-3)
{ return E*E/wigner_seitz_radius(n); } 

// debye length of a particular plasma species, cgs units
double debye_length(double n,double T) 
// n: plasma density (cm^-3)
// T: temperature (K) of a particular plasma species
{ return sqrt(K_B*T/(4*PI*n*E*E));}

// plasma screening parameter for electrons, dimensionless, inputs cgs
double screening_parameter(double n, double Te)
// n: plasma density (cm^-3)
// Te: electron temperature (K)
{ return wigner_seitz_radius(n)/debye_length(n,Te); } 

// equilibrium plasma temperature due to disorder-induced heating, inputs cgs
double dih_temp(double n,double Te) 
// n: plasma density (cm^-3)
// Te: electron temperature (K)
{ return 2.*nearest_coulomb_pot(n)/3./K_B*(1.+screening_parameter(n,Te)/2.); }

// timescale for a Gaussian UCNP expansion into vacuum, cgs units
double tau_exp(double sig,double m, double Te)
// n: plasma density (cm^-3)
// m: ion mass (g)
// Te: electron temperature (K)
{
    return sqrt(m*pow(sig,2.)/(K_B*Te));
}
}