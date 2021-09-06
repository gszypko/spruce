#include "solarutils.hpp"

namespace SolarUtils {

  MhdInp SolarMHDInput(const PlasmaSettings& pms)
  {
    std::string problem_type = pms.getopt("type");
    const int xdim = (int)(pms.getvar("xdim")+0.01);
    const int ydim = (int)(pms.getvar("ydim")+0.01);
    const double init_temp = pms.getvar("init_temp");
    // const double x_size = pms.getvar("x_size");
    // const double y_size = pms.getvar("y_size");
    const double b_0 = pms.getvar("b_0");
    const double n_base = pms.getvar("n_base");
    double ion_mass = pms.getvar("ion_mass");
    double adiabatic_index = pms.getvar("adiabatic_index");
    double duration = pms.getvar("duration");

    double base_rho = ion_mass*n_base; //initial mass density at base in g cm ^-3
    double scale_height = 2.0*K_B*init_temp/(ion_mass*BASE_GRAV);
    Grid pos_x(xdim,ydim); Grid pos_y(xdim,ydim);
    double dx = 4.0*PI*scale_height/xdim, dy = 4.0*PI*scale_height/ydim; //this ensures correct periodicity of bipolar field
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        pos_x(i,j) = (i - (double)(xdim-1)*0.5)*dx;
        pos_y(i,j) = j*dy;
      }
    }

    MhdInp mi(xdim,ydim);
    mi.set_var(PlasmaDomain::pos_x,pos_x);
    mi.set_var(PlasmaDomain::pos_y,pos_y);
    mi.set_var(PlasmaDomain::mom_x, Grid(xdim,ydim,0.0));
    mi.set_var(PlasmaDomain::mom_y, Grid(xdim,ydim,0.0));

    if(problem_type == "corona"){
      mi.set_var(PlasmaDomain::temp, Grid(xdim,ydim,init_temp));
      mi.set_var(PlasmaDomain::rho, HydrostaticFalloff(base_rho,scale_height,pos_y));
      mi.set_var(PlasmaDomain::b_x, BipolarField(pos_x, pos_y, b_0, scale_height, 0));
      mi.set_var(PlasmaDomain::b_y, BipolarField(pos_x, pos_y, b_0, scale_height, 1));
      mi.set_var(PlasmaDomain::b_z, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::grav_y, SolarGravity(BASE_GRAV,R_SUN,pos_y));
      mi.set_var(PlasmaDomain::grav_x, Grid(xdim,ydim,0.0));
    } else if(problem_type == "uniform"){
      mi.set_var(PlasmaDomain::temp, Grid(xdim,ydim,init_temp));
      mi.set_var(PlasmaDomain::rho, Grid(xdim,ydim,base_rho*std::exp(-2.0*PI)));
      mi.set_var(PlasmaDomain::b_x, BipolarField(pos_x, pos_y, b_0, scale_height, 0));
      mi.set_var(PlasmaDomain::b_y, BipolarField(pos_x, pos_y, b_0, scale_height, 1));
      mi.set_var(PlasmaDomain::b_z, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::grav_y, SolarGravity(BASE_GRAV,R_SUN,pos_y));
      mi.set_var(PlasmaDomain::grav_x, Grid(xdim,ydim,0.0));
    } else if(problem_type == "gaussian"){
      double stdevx = pms.getvar("stdevx");
      double stdevy = pms.getvar("stdevy");
      // mi.set_var(PlasmaDomain::temp, GaussianGrid(xdim,ydim,init_temp,10.0*init_temp,xdim/stdevx,ydim/stdevy));
      mi.set_var(PlasmaDomain::temp, Grid(xdim,ydim,init_temp));
      mi.set_var(PlasmaDomain::rho, GaussianGrid(xdim,ydim,0.1*base_rho,10.0*base_rho,xdim/stdevx,ydim/stdevy));
      mi.set_var(PlasmaDomain::b_x, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::b_y, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::b_z, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::grav_y, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::grav_x, Grid(xdim,ydim,0.0));
    }

    mi.set_ion_mass(ion_mass);
    mi.set_adiabatic_index(adiabatic_index);
    mi.set_duration(duration);

    return mi;
  }


  //Sets gravity to fall off from base_gravity at bottom of the domain,
  //as though from the surface of a planet with radius r_solar
  Grid SolarGravity(double base_gravity, double r_solar, const Grid &m_pos_y)
  {
    #if BENCHMARKING_ON
    InstrumentationTimer timer(__PRETTY_FUNCTION__);
    #endif
    int m_xdim = m_pos_y.rows(), m_ydim = m_pos_y.cols();
    Grid result(m_xdim,m_ydim);
    for(int j=0; j<m_ydim; j++){
      for(int i=0; i<m_xdim; i++){
        result(i,j) = -base_gravity*std::pow(r_solar/(r_solar+m_pos_y(i,j)),2.0);
      }
    }
    return result;
  }

  //Generates potential bipolar field for component corresponding to index "index"
  //Centered s.t. origin lies at bottom middle of domain
  //Pressure scale height h, field poles at +/- l, field strength at poles b0
  Grid BipolarField(const Grid& m_pos_x, const Grid& m_pos_y, double b0, double h, int index){
    #if BENCHMARKING_ON
    InstrumentationTimer timer(__PRETTY_FUNCTION__);
    #endif
    int m_xdim = m_pos_x.rows(); int m_ydim = m_pos_y.cols();
    assert(m_xdim == m_ydim && "pos_x and pos_y must have same dimensions");
    Grid result = Grid::Zero(m_xdim, m_ydim);
    #pragma omp parallel for collapse(2)
    for(int i=0; i<m_xdim; i++){
      for(int j=0; j<m_ydim; j++){
        if(index == 0) result(i,j) = b0*std::exp(-0.5*m_pos_y(i,j)/h)*std::cos(0.5*m_pos_x(i,j)/h);
        else result(i,j) = -b0*std::exp(-0.5*m_pos_y(i,j)/h)*std::sin(0.5*m_pos_x(i,j)/h);
      }
    }
    return result;
  }

  //Generates grid with hydrostatic falloff in the y-direction, with the quantity
  //"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
  //Assumes that gravity is set to vary with distance.
  Grid HydrostaticFalloff(double base_value, double scale_height, const Grid& m_pos_y){
    #if BENCHMARKING_ON
    InstrumentationTimer timer(__PRETTY_FUNCTION__);
    #endif
    Grid result = Grid::Zero(m_pos_y.rows(), m_pos_y.cols());
    int xdim = m_pos_y.rows(), ydim = m_pos_y.cols();
    #pragma omp parallel for collapse(2)
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        //result(i,j) = base_value*std::exp(-y/scale_height);
        result(i,j) = base_value*std::exp(M_SUN*GRAV_CONST/BASE_GRAV/scale_height*(1.0/(m_pos_y(i,j)+R_SUN) - 1.0/R_SUN));
      }
    }
    return result;
  }

  //Generates gaussian initial condition for a variable, centered at middle of grid
  //std_dev_x and std_dev_y are the standard deviation of the distribution in the x
  //and y directions, in units of grid cell widths
  Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y){
    #if BENCHMARKING_ON
    InstrumentationTimer timer(__PRETTY_FUNCTION__);
    #endif
    std::vector<double> gauss_x(xdim), gauss_y(ydim);
    for(int i=0; i<xdim; i++){
      gauss_x[i] = std::exp(-0.5*std::pow(((double)i-0.5*(double)(xdim-1))/std_dev_x,2.0));
    }
    for(int j=0; j<ydim; j++){
      gauss_y[j] = std::exp(-0.5*std::pow(((double)j-0.5*(double)(ydim-1))/std_dev_y,2.0));
    }
    Grid result(xdim,ydim);
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        result(i,j) = (max-min)*gauss_x[i]*gauss_y[j] + min;
      }
    }
    return result;
  }

}
