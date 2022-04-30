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
    const double n_base = pms.getvar("n_base");
    double ion_mass = pms.getvar("ion_mass");
    double adiabatic_index = pms.getvar("adiabatic_index");
    double duration = pms.getvar("duration");
    double growth_factor = pms.getvar("growth_factor");

    double base_rho = ion_mass*n_base; //initial mass density at base in g cm ^-3
    double scale_height = 2.0*K_B*init_temp/(ion_mass*BASE_GRAV);
    Grid pos_x(xdim,ydim); Grid pos_y(xdim,ydim);
    Grid d_x(xdim,ydim); Grid d_y(xdim,ydim);
    double dx = 4.0*PI*scale_height/xdim, dy = 4.0*PI*scale_height/ydim; //this ensures correct periodicity of bipolar field

    if(problem_type == "corona" || problem_type == "uniform"){
      dy = scale_height*pms.getvar("y_size")/ydim;
      double uniform_fraction = pms.getvar("uniform_fraction");
      double curr_scaling = 1.0;
      for(int i=0; i<xdim; i++)
        for(int j=0; j<ydim; j++) d_x(i,j) = dx;
      for(int j=0; j<ydim; j++){
        if(j >= uniform_fraction*ydim) curr_scaling *= growth_factor;
        for(int i=0; i<xdim; i++) d_y(i,j) = curr_scaling*dy;
      }
    } else {
      assert(xdim%4 == 0 && ydim%4 == 0);
      int x_core_size = xdim/2, y_core_size = ydim/2;
      for(int j=0; j<ydim; j++){
        for(int i=(xdim/2 - x_core_size/2); i<(xdim/2 + x_core_size/2); i++){
          d_x(i,j) = dx;
        }
        double this_dx = growth_factor*dx;
        for(int i=(xdim/2 - x_core_size/2) - 1; i>=0; i--){
          d_x(i,j) = this_dx;
          this_dx *= growth_factor;
        }
        this_dx = growth_factor*dx;
        for(int i=(xdim/2 + x_core_size/2); i<xdim; i++){
          d_x(i,j) = this_dx;
          this_dx *= growth_factor;
        }
      }
      for(int i=0; i<xdim; i++){
        for(int j=(ydim/2 - y_core_size/2); j<(ydim/2 + y_core_size/2); j++){
          d_y(i,j) = dy;
        }
        double this_dy = growth_factor*dy;
        for(int j=(ydim/2 - y_core_size/2) - 1; j>=0; j--){
          d_y(i,j) = this_dy;
          this_dy *= growth_factor;
        }
        this_dy = growth_factor*dy;
        for(int j=(ydim/2 + y_core_size/2); j<ydim; j++){
          d_y(i,j) = this_dy;
          this_dy *= growth_factor;
        }
      }
    }

    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        if(i==0) pos_x(i,j) = 0.5*d_x(i,j);
        else pos_x(i,j) = pos_x(i-1,j) + 0.5*d_x(i-1,j) + 0.5*d_x(i,j);
        if(j==0) pos_y(i,j) = 0.5*d_y(i,j);
        else pos_y(i,j) = pos_y(i,j-1) + 0.5*d_y(i,j-1) + 0.5*d_y(i,j);
      }
    }

    std::cout << "woo\n";
    MhdInp mi(xdim,ydim);
    mi.set_var(PlasmaDomain::pos_x, pos_x);
    mi.set_var(PlasmaDomain::pos_y, pos_y);
    mi.set_var(PlasmaDomain::mom_x, Grid(xdim,ydim,0.0));
    mi.set_var(PlasmaDomain::mom_y, Grid(xdim,ydim,0.0));
    mi.set_var(PlasmaDomain::bi_x, Grid(xdim,ydim,0.0));
    mi.set_var(PlasmaDomain::bi_y, Grid(xdim,ydim,0.0));
    std::cout << "foo\n";

    if(problem_type == "corona"){
      double b_0 = pms.getvar("b_0");
      mi.set_var(PlasmaDomain::d_x,d_x);
      mi.set_var(PlasmaDomain::d_y,d_y);
      mi.set_var(PlasmaDomain::temp, Grid(xdim,ydim,init_temp));
      Grid b_x, b_y;
      b_x = BipolarField(pos_x, pos_y, b_0, scale_height, 0);
      b_y = BipolarField(pos_x, pos_y, b_0, scale_height, 1);
      std::cout << "loo\n";
      mi.set_var(PlasmaDomain::be_x, b_x);
      mi.set_var(PlasmaDomain::be_y, b_y);
      mi.set_var(PlasmaDomain::grav_y, SolarGravity(BASE_GRAV,R_SUN,pos_y));
      mi.set_var(PlasmaDomain::grav_x, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::rho, HydrostaticFalloff(base_rho,scale_height,pos_y));
    } else if(problem_type == "gaussian"){
      mi.set_var(PlasmaDomain::d_x,d_x);
      mi.set_var(PlasmaDomain::d_y,d_y);
      double stdevx = pms.getvar("stdevx");
      double stdevy = pms.getvar("stdevy");
      double b_x = pms.getvar("b_x");
      double b_y = pms.getvar("b_y");
      mi.set_var(PlasmaDomain::temp, Grid(xdim,ydim,init_temp));
      mi.set_var(PlasmaDomain::rho, GaussianGrid(pos_x,pos_y,0.1*base_rho,10.0*base_rho,xdim/stdevx,ydim/stdevy));
      mi.set_var(PlasmaDomain::be_x, Grid(xdim,ydim,b_x));
      mi.set_var(PlasmaDomain::be_y, Grid(xdim,ydim,b_y));
      mi.set_var(PlasmaDomain::grav_y, Grid(xdim,ydim,0.0));
      mi.set_var(PlasmaDomain::grav_x, Grid(xdim,ydim,0.0));
    }

    mi.set_ion_mass(ion_mass);
    mi.set_adiabatic_index(adiabatic_index);
    mi.set_duration(duration);

    std::cout << "noo\n";
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
        result(i,j) = -base_gravity;
        //result(i,j) = -base_gravity*std::pow(r_solar/(r_solar+m_pos_y(i,j)),2.0);
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
    // assert(m_xdim == m_ydim && "pos_x and pos_y must have same dimensions");
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
        result(i,j) = base_value*std::exp(-m_pos_y(i,j)/scale_height);
        //result(i,j) = base_value*std::exp(M_SUN*GRAV_CONST/BASE_GRAV/scale_height*(1.0/(m_pos_y(i,j)+R_SUN) - 1.0/R_SUN));
      }
    }
    return result;
  }

  //Generates gaussian initial condition for a variable, centered at middle of grid
  //std_dev_x and std_dev_y are the standard deviation of the distribution in the x
  //and y directions (units of distance)
  Grid GaussianGrid(const Grid& pos_x, const Grid& pos_y, double min, double max, double std_dev_x, double std_dev_y){
    assert(pos_x.rows() == pos_y.rows() && pos_x.cols() == pos_y.cols());
    int xdim = pos_x.rows(), ydim = pos_x.cols();
    return GaussianGrid(xdim, ydim, min, max, std_dev_x, std_dev_y, 0.5*(double)(xdim-1), 0.5*(double)(ydim-1));
  }

  Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y,
                    double center_x, double center_y){
    #if BENCHMARKING_ON
    InstrumentationTimer timer(__PRETTY_FUNCTION__);
    #endif
    std::vector<double> gauss_x(xdim), gauss_y(ydim);
    for(int i=0; i<xdim; i++){
      gauss_x[i] = std::exp(-0.5*std::pow(((double)i-center_x)/std_dev_x,2.0));
    }
    for(int j=0; j<ydim; j++){
      gauss_y[j] = std::exp(-0.5*std::pow(((double)j-center_y)/std_dev_y,2.0));
    }
    Grid result(xdim,ydim);
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        result(i,j) = (max-min)*gauss_x[i]*gauss_y[j] + min;
      }
    }
    // std::cout << result << std::endl;
    return result;
  }


}
