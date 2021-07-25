#include "constants.hpp"
#include "utils.hpp"
#include "derivs.hpp"
#include "grid.hpp"
#include <omp.h>
#include <cmath>
#include <vector>
#include <limits>

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

//Generates gaussian initial condition for a variable, centered at middle of grid
//std_dev_x and std_dev_y are the standard deviation of the distribution in the x
//and y directions, in units of grid cell widths
Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  std::vector<double> gauss_x(xdim), gauss_y(ydim);
  // double sigmax = 0.05*xdim;
  // double sigmay = 0.05*ydim;
  for(int i=0; i<xdim; i++){
    gauss_x[i] = std::exp(-0.5*std::pow(((double)i-0.5*(double)(xdim-1))/std_dev_x,2.0));
  }
  for(int j=0; j<ydim; j++){
    gauss_y[j] = std::exp(-0.5*std::pow(((double)j-0.5*(double)(ydim-1))/std_dev_y,2.0));
  }
  Grid result(xdim,ydim);
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      result(i,j) = gauss_x[i]*gauss_y[j];
    }
  }
  return result;
}

//Generates potential bipolar field for component corresponding to index "index"
//Centered s.t. origin lies at bottom middle of domain
//Pressure scale height h, field poles at +/- l, field strength at poles b0
Grid BipolarField(int xdim, int ydim, double b0, double h, double dx, double dy, int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result = Grid::Zero(xdim, ydim);
  #pragma omp parallel for collapse(2)
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      double x = (i - (double)(xdim-1)*0.5)*dx;
      double y = j*dy;
      if(index == 0) result(i,j) = b0*std::exp(-0.5*y/h)*std::cos(0.5*x/h);
      else result(i,j) = -b0*std::exp(-0.5*y/h)*std::sin(0.5*x/h);
    }
  }
  return result;
}

//Generates grid with exponential falloff in the y-direction, with the quantity
//"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
Grid HydrostaticFalloff(double base_value, double scale_height, int xdim, int ydim, double dy){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result = Grid::Zero(xdim, ydim);
  #pragma omp parallel for collapse(2)
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      double y = j*dy;
      //result(i,j) = base_value*std::exp(-y/scale_height);
      result(i,j) = base_value*std::exp(M_SUN*GRAV_CONST/BASE_GRAV/scale_height*(1/(y+R_SUN) - 1/R_SUN));
    }
  }
  return result;
}

//Computes 1D cell-centered conductive flux from temperature "temp"
//Flux computed in direction indicated by "index": 0 for x, 1 for y
//k0 is conductive coefficient
Grid onedim_conductive_flux(const Grid &temp, const Grid &rho, double k0, int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid kappa_max = DX*DY*K_B*(rho/M_I)/DT_THERMAL_MIN;
  int xdim = temp.rows();
  int ydim = temp.cols();
  Grid flux = Grid::Zero(xdim,ydim);
  flux = temp.pow(7.0/2.0);
  #pragma omp parallel for collapse(2)
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      flux(i,j) = std::pow(temp(i,j),7.0/2.0);
    }
  }
  return -(k0*temp.pow(5.0/2.0)).min(kappa_max)*derivative1D(temp,index);
  // return -2.0/7.0*k0*derivative1D(flux,index);
}

//Computes cell-centered, field-aligned conductive flux from temperature "temp"
//temp is temperature Grid
//b_hat_x, b_hat_y are the components of the *unit* vector b_hat
//k0 is conductive coefficient
//Output is written to flux_out_x and flux_out_y
void field_aligned_conductive_flux(Grid &flux_out_x, Grid &flux_out_y, const Grid &temp, const Grid &rho,
                                    const Grid &b_hat_x, const Grid &b_hat_y, double k0){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  int xdim = temp.rows();
  int ydim = temp.cols();
  Grid con_flux_x = onedim_conductive_flux(temp, rho, k0, 0);
  Grid con_flux_y = onedim_conductive_flux(temp, rho, k0, 1);
  #pragma omp parallel
  {
    #if BENCHMARKING_ON
    InstrumentationTimer timer("loop thread");
    #endif
    #pragma omp for collapse(2)
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        double flux_magnitude = con_flux_x(i,j)*b_hat_x(i,j) + con_flux_y(i,j)*b_hat_y(i,j);
        flux_out_x(i,j) = flux_magnitude*b_hat_x(i,j);
        flux_out_y(i,j) = flux_magnitude*b_hat_y(i,j);
      }
    }
  }
}

//Computes saturated conductive flux at each point in grid,
//then ensures that provided fluxes do not exceed the saturation point
void saturate_conductive_flux(Grid &flux_out_x, Grid &flux_out_y, const Grid &rho, const Grid &temp){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid sat_flux_mag = (1.0/6.0)*(3.0/2.0)*(rho/M_I)*(K_B*temp).pow(1.5)/std::sqrt(M_ELECTRON);
  Grid flux_mag = (flux_out_x.square() + flux_out_y.square()).sqrt();
  Grid scale_factor = sat_flux_mag /((sat_flux_mag.square() + flux_mag.square()).sqrt());
  flux_out_x *= scale_factor;
  flux_out_y *= scale_factor;
}
