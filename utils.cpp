#include "constants.hpp"
#include "utils.hpp"
#include "derivs.hpp"
#include "grid.hpp"
#include <omp.h>
#include <cmath>
#include <limits>

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

//Enforce dynamic time stepping
//Local dt values returned in local_dt
// double recompute_dt(const Grid &press, const Grid &rho, const Grid &vx, const Grid &vy, Grid &local_dt){
//   #if BENCHMARKING_ON
//   InstrumentationTimer timer(__PRETTY_FUNCTION__);
//   #endif
//   Grid c_s = (GAMMA*press/rho).sqrt();
//   Grid abs_vx = DX/(c_s+vx.abs());
//   Grid abs_vy = DY/(c_s+vy.abs());
//   local_dt = abs_vx.min(abs_vy);
//   return local_dt.min();
// }

//Enforce dynamic time stepping for thermal conduction
//Local dt values returned in local_dt
// double recompute_dt_thermal(const Grid &rho, const Grid &temp, Grid &local_dt, const Grid &b_hat_x, const Grid &b_hat_y){
//   #if BENCHMARKING_ON
//   InstrumentationTimer timer(__PRETTY_FUNCTION__);
//   #endif
//   if(DT_CALCULATION_MODE == 0) local_dt = K_B/KAPPA_0*(rho/M_I)*DX*DY/temp.pow(2.5);
//   else{
//     Grid kappa_modified;
//     Grid field_temp_gradient = derivative1D(temp,0)*b_hat_x + derivative1D(temp,1)*b_hat_y;
//     Grid con_flux_x = Grid::Zero(XDIM,YDIM);
//     Grid con_flux_y = Grid::Zero(XDIM,YDIM);
//     field_aligned_conductive_flux(con_flux_x, con_flux_y, temp, rho, b_hat_x, b_hat_y, KAPPA_0);
//     Grid flux_c = (con_flux_x.square() + con_flux_y.square()).sqrt();
//     saturate_conductive_flux(con_flux_x, con_flux_y, rho, temp);
//     Grid flux_sat = (con_flux_x.square() + con_flux_y.square()).sqrt();
//     if(DT_CALCULATION_MODE == 1){ //coefficent saturation
//       kappa_modified = (flux_sat/temp.pow(2.5)/field_temp_gradient).abs();
//       local_dt = K_B/kappa_modified*(rho/M_I)*DX*DY/temp.pow(2.5);
//     }
//     else if(DT_CALCULATION_MODE == 2){ //whole-term saturation
//       kappa_modified = (flux_sat/field_temp_gradient).abs();
//       local_dt = K_B/kappa_modified*(rho/M_I)*DX*DY;
//     } 
//   }
//   return local_dt.min();
// }

//Enforce minimum dynamic time step for radiation
//Local dt values returned in local_dt (negative values indicate no radiation locally occuring)
// double recompute_dt_radiative(const Grid &energy, const Grid &rad_loss_rate, Grid &local_dt){
//   #if BENCHMARKING_ON
//   InstrumentationTimer timer(__PRETTY_FUNCTION__);
//   #endif
//   local_dt = Grid::Zero(XDIM,YDIM);
//   double running_min_dt = std::numeric_limits<double>::max();
//   #pragma omp parallel for reduction(min: running_min_dt) collapse(2)
//   for(int i=0; i<XDIM; i++){
//     for(int j=0; j<YDIM; j++){
//       if(rad_loss_rate(i,j) > 0.0){
//         double this_dt = std::abs(energy(i,j)/rad_loss_rate(i,j));
//         local_dt(i,j) = this_dt;
//         running_min_dt = std::min(running_min_dt, this_dt);
//       }
//     }
//   }
//   return running_min_dt;
// }

//Generates gaussian initial condition for a variable, centered at middle of grid
// Grid GaussianGrid(int xdim, int ydim, double min, double max){
//   #if BENCHMARKING_ON
//   InstrumentationTimer timer(__PRETTY_FUNCTION__);
//   #endif
//   Eigen::VectorXd gauss_x(xdim), gauss_y(ydim);
//   double sigmax = 0.05*xdim;
//   double sigmay = 0.05*ydim;
//   for(int i=0; i<xdim; i++){
//     gauss_x(i) = std::exp(-0.5*std::pow(((double)i-0.5*(double)(xdim-1))/sigmax,2.0));
//   }
//   for(int j=0; j<ydim; j++){
//     gauss_y(j) = std::exp(-0.5*std::pow(((double)j-0.5*(double)(ydim-1))/sigmay,2.0));
//   }
//   Eigen::MatrixXd gauss_2d = gauss_x * gauss_y.transpose();
//   gauss_2d = gauss_2d*(max-min) + min*Eigen::MatrixXd::Ones(xdim,ydim);
//   return gauss_2d.array();
// }

//Generates potential bipolar field for component corresponding to index "index"
//Centered s.t. origin lies at bottom middle of domain
//Pressure scale height h, field poles at +/- l, field strength at poles b0
Grid BipolarField(const int xdim, const int ydim, const double b0, const double h, const int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result = Grid::Zero(xdim, ydim);
  #pragma omp parallel for collapse(2)
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      double x = (i - (double)(xdim-1)*0.5)*DX;
      double y = j*DY;
      if(index == 0) result(i,j) = b0*std::exp(-0.5*y/h)*std::cos(0.5*x/h);
      else result(i,j) = -b0*std::exp(-0.5*y/h)*std::sin(0.5*x/h);
    }
  }
  return result;
}

//Generates grid with exponential falloff in the y-direction, with the quantity
//"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
Grid HydrostaticFalloff(const double base_value, const double scale_height, const int xdim, const int ydim){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result = Grid::Zero(xdim, ydim);
  #pragma omp parallel for collapse(2)
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      double y = j*DY;
      result(i,j) = base_value*std::exp(-y/scale_height);
    }
  }
  return result;
}

//Generates Grid containing magnitude of gravitational
//acceleration (in y-direction) at each grid cell
// Grid Gravity(const double base_grav, const double r_sun, const int xdim, const int ydim){
//   #if BENCHMARKING_ON
//   InstrumentationTimer timer(__PRETTY_FUNCTION__);
//   #endif
//   Grid result = Grid::Zero(xdim,ydim);
//   for(int j=0; j<ydim; j++){
//     double y = j*DY;
//     for(int i=0; i<xdim; i++){
//       result(i,j) = base_grav*std::pow(r_sun/(r_sun+y),2.0);
//     }
//   }
//   return result;
// }

//Computes 1D cell-centered conductive flux from temperature "temp"
//Flux computed in direction indicated by "index": 0 for x, 1 for y
//k0 is conductive coefficient
Grid onedim_conductive_flux(const Grid &temp, const Grid &rho, const double k0, const int index){
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
                                    const Grid &b_hat_x, const Grid &b_hat_y, const double k0){
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

// Grid radiative_losses(const Grid &rho, const Grid &temp, const int xdim, const int ydim){
//   #if BENCHMARKING_ON
//   InstrumentationTimer timer(__PRETTY_FUNCTION__);
//   #endif
//   Grid result = Grid::Zero(xdim,ydim);
//   #pragma omp parallel
//   {
//     #if BENCHMARKING_ON
//     InstrumentationTimer timer("loop thread");
//     #endif
//     #pragma omp for collapse(2)
//     for(int i=0; i<xdim; i++){
//       for(int j=0; j<ydim; j++){
//         if(temp(i,j) < TEMP_CHROMOSPHERE){
//           result(i,j) = 0.0;
//           continue;
//         }
//         double logtemp = std::log10(temp(i,j));
//         double n = rho(i,j)/M_I;
//         double chi, alpha;
//         if(logtemp <= 4.97){
//           chi = 1.09e-31;
//           alpha = 2.0;
//         } else if(logtemp <= 5.67){
//           chi = 8.87e-17;
//           alpha = -1.0;
//         } else if(logtemp <= 6.18){
//           chi = 1.90e-22;
//           alpha = 0.0;
//         } else if(logtemp <= 6.55){
//           chi = 3.53e-13;
//           alpha = -1.5;
//         } else if(logtemp <= 6.90){
//           chi = 3.46e-25;
//           alpha = 1.0/3.0;
//         } else if(logtemp <= 7.63){
//           chi = 5.49e-16;
//           alpha = -1.0;
//         } else{
//           chi = 1.96e-27;
//           alpha = 0.5;
//         }
//         result(i,j) = n*n*chi*std::pow(temp(i,j),alpha);
//         if(temp(i,j) < TEMP_CHROMOSPHERE + RADIATION_RAMP){
//           //LINEAR RAMP
//           // double ramp = (temp(i,j) - TEMP_CHROMOSPHERE)/RADIATION_RAMP;
//           //COSINE RAMP
//           double ramp = 0.5*(1.0 - std::cos((temp(i,j) - TEMP_CHROMOSPHERE)*PI/RADIATION_RAMP));
//           result(i,j) *= ramp;
//         }
//       }
//     }
//   }
//   return result;
// }

