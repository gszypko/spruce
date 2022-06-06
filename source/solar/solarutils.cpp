#include "solarutils.hpp"
#include "constants.hpp"
#include "plasmadomain.hpp"
#include <cmath>

namespace SolarUtils {

  //Sets gravity to fall off from base_gravity at bottom of the domain,
  //as though from the surface of a planet with radius r_solar
  Grid SolarGravity(double base_gravity, double r_solar, const Grid &m_pos_y)
  {
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
