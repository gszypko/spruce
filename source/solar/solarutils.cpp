#include <cmath>
#include "grid.hpp"
#include "solarutils.hpp"
#include "constants.hpp"

namespace SolarUtils {
  static const double temp_chromosphere = 3.0e4;
  static const int xdim = 100;
  static const int ydim = 100;
  static const double dx = 2.2649e7;
  static const double dy = 2.2649e7;
  static const double b_0 = 100.0;

  void SolarInitialize(Grid& rho, Grid& temp, Grid& mom_x, Grid& mom_y,
                       Grid& b_x, Grid& b_y, Grid& b_z,
                       Grid& pos_x, Grid& pos_y, Grid& grav_x, Grid& grav_y)
  {
    rho = Grid(xdim,ydim,0.0);
    temp = Grid(xdim,ydim,0.0);
    mom_x = Grid(xdim,ydim,0.0); mom_y = Grid(xdim,ydim,0.0);
    b_x = Grid(xdim,ydim,0.0); b_y = Grid(xdim,ydim,0.0); b_z = Grid(xdim,ydim,0.0);
    pos_x = Grid(xdim,ydim,0.0); pos_y = Grid(xdim,ydim,0.0);
    grav_x = Grid(xdim,ydim,0.0); grav_y = Grid(xdim,ydim,0.0);
    double base_rho = M_I*1.0e12; //initial mass density at base, g cm^-3
    double scale_height = 2.0*K_B*temp_chromosphere/(M_I*BASE_GRAV);
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        pos_x(i,j) = (i - (double)(xdim-1)*0.5)*dx;
        pos_y(i,j) = j*dy;
      }
    }
    rho = HydrostaticFalloff(base_rho,scale_height,pos_y);
    temp = Grid(xdim,ydim,temp_chromosphere);
    b_x = BipolarField(pos_x, pos_y, b_0, scale_height, 0);
    b_y = BipolarField(pos_x, pos_y, b_0, scale_height, 1);
    grav_y = SolarGravity(BASE_GRAV,R_SUN,pos_y);
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
    #pragma omp parallel for collapse(2)
    for(int i=0; i<m_pos_y.rows(); i++){
      for(int j=0; j<m_pos_y.cols(); j++){
        //result(i,j) = base_value*std::exp(-y/scale_height);
        result(i,j) = base_value*std::exp(M_SUN*GRAV_CONST/BASE_GRAV/scale_height*(1.0/(m_pos_y(i,j)+R_SUN) - 1.0/R_SUN));
      }
    }
    return result;
  }
}