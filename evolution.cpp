//evolution.cpp
//PlasmaDomain functionality relating to time-evolution of the plasma

#include <cmath>
#include <omp.h>
#include "utils.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"

void PlasmaDomain::run()
{
  outputPreamble();
  outputCurrentState();
  while( (max_iterations < 0 || m_iter < max_iterations) && (max_time < 0.0 || m_time < max_time) ){
    #if BENCHMARKING_ON
    InstrumentationTimer timer((std::string("iteration ") + std::to_string(m_iter)).c_str());
    #endif
    double old_time = m_time;
    advanceTime();
    if((iter_output_interval > 0 && m_iter%iter_output_interval == 0) || (time_output_interval > 0.0 &&
       (int)(m_time/time_output_interval) > (int)(old_time/time_output_interval))) outputCurrentState();
    writeStateFile();
  }
  cleanUpStateFiles();
}

void PlasmaDomain::advanceTime(bool verbose)
{
  double min_dt, min_dt_thermal, min_dt_rad;
  int subcycles_thermal, subcycles_rad;

  double dt_raw = m_grids[dt].min();
  min_dt = epsilon*dt_raw;

  if(thermal_conduction){
    min_dt_thermal = std::max(epsilon_thermal*m_grids[dt_thermal].min(),dt_thermal_min);
    subcycles_thermal = (int)(min_dt/min_dt_thermal)+1;
  }
  if(radiative_losses){
    min_dt_rad = epsilon_rad*m_grids[dt_rad].min();
    subcycles_rad = (int)(min_dt/min_dt_rad)+1;
  }

  if(verbose) printUpdate(min_dt,subcycles_thermal,subcycles_rad);

  if(thermal_conduction) subcycleConduction(subcycles_thermal,min_dt);
  if(radiative_losses) subcycleRadiation(subcycles_rad,min_dt);

  //Advance time by min_dt
  Grid &m_mom_x = m_grids[mom_x], &m_mom_y = m_grids[mom_y],
        &m_v_x = m_grids[v_x], &m_v_y = m_grids[v_y],
        &m_rho = m_grids[rho], &m_energy = m_grids[energy], &m_press = m_grids[press];
  Grid viscous_force_x = epsilon_viscous*(0.5*m_dx*m_dx/dt_raw)*laplacian(m_mom_x);
  Grid viscous_force_y = epsilon_viscous*(0.5*m_dy*m_dy/dt_raw)*laplacian(m_mom_y);
  Grid zero(1,1,0.0);
  Grid rho_next = m_rho - min_dt*divergence(m_rho,zero,zero);
  Grid mom_x_next = m_mom_x - min_dt*divergence(m_mom_x, m_press + m_grids[mag_pxx], m_grids[mag_pxy])
                    + min_dt*m_rho*m_grids[grav_x] + min_dt*viscous_force_x;
  Grid mom_y_next = m_mom_y - min_dt*divergence(m_mom_y, m_grids[mag_pxy], m_press + m_grids[mag_pyy])
                    + min_dt*m_rho*m_grids[grav_y] + min_dt*viscous_force_y;
  Grid energy_next = m_energy - min_dt*divergence(m_energy+m_press, m_grids[mag_pxx]*m_v_x + m_grids[mag_pxy]*m_v_y, m_grids[mag_pxy]*m_v_x + m_grids[mag_pyy]*m_v_y) 
                    + min_dt*m_rho*(m_v_x*m_grids[grav_x] + m_v_y*m_grids[grav_y]) + min_dt*(m_v_x*viscous_force_x + m_v_y*viscous_force_y);
  if(ambient_heating) energy_next += min_dt*heating_rate;

  clampWallBoundaries(mom_x_next, mom_y_next, rho_next, energy_next);

  m_rho = rho_next;
  m_mom_x = mom_x_next;
  m_mom_y = mom_y_next;
  m_energy = energy_next;
  
  m_time += min_dt;
  m_iter++;

  catchUnderdensity();
  recomputeTemperature();
  recomputeDerivedVariables();
}

void PlasmaDomain::subcycleConduction(int subcycles_thermal, double dt_total)
{
  //Subcycle to simulate field-aligned thermal conduction
  // Grid energy_relaxed = energy;
  Grid &m_energy = m_grids[energy], &m_temp = m_grids[temp], &m_press = m_grids[press];
  Grid nonthermal_energy = 0.5*(m_grids[mom_x]*m_grids[v_x] + m_grids[mom_y]*m_grids[v_y]) + m_grids[mag_press];
  Grid energy_floor = nonthermal_energy + thermal_energy_min; //To ensure non-negative thermal pressure
  Grid energy_next;
  Grid dummy_grid(m_xdim,m_ydim,0.0); //to feed into clampWallBoundaries, since rho and mom are unchanged here
  for(int subcycle = 0; subcycle < subcycles_thermal; subcycle++){
    Grid con_flux_x(m_xdim,m_ydim,0.0);
    Grid con_flux_y(m_xdim,m_ydim,0.0);
    fieldAlignedConductiveFlux(con_flux_x, con_flux_y, m_temp, m_grids[rho], m_grids[b_hat_x], m_grids[b_hat_y], KAPPA_0);
    if(flux_saturation) saturateConductiveFlux(con_flux_x, con_flux_y, m_grids[rho], m_temp);
    energy_next = m_energy - (dt_total/(double)subcycles_thermal)*(derivative1D(con_flux_x,0)+derivative1D(con_flux_y,1));
    clampWallBoundaries(dummy_grid, dummy_grid, dummy_grid, energy_next);
    m_energy = energy_next.max(energy_floor);
    m_press = (GAMMA - 1.0)*(m_energy - nonthermal_energy);
    m_temp = M_I*m_press/(2.0*K_B*m_grids[rho]);
    // Enforce temp floor everywhere
    m_temp = m_temp.max(temp_min);
    // Recompute energy after flooring temperature
    m_press = 2.0*K_B*m_grids[rho]*m_temp/M_I;
    m_energy = m_press/(GAMMA - 1.0) + nonthermal_energy;
  }
}

void PlasmaDomain::subcycleRadiation(int subcycles_rad, double dt_total)
{
  Grid &m_energy = m_grids[energy], &m_temp = m_grids[temp], &m_press = m_grids[press];
  //Subcycle to simulate radiative losses
  Grid nonthermal_energy = 0.5*(m_grids[mom_x]*m_grids[v_x] + m_grids[mom_y]*m_grids[v_y]) + m_grids[mag_press];
  Grid energy_floor = nonthermal_energy + thermal_energy_min; //To ensure non-negative thermal pressure
  Grid energy_next;
  Grid dummy_grid(m_xdim,m_ydim,0.0); //to feed into clampWallBoundaries, since rho and mom are unchanged here
  for(int subcycle = 0; subcycle < subcycles_rad; subcycle++){
    recomputeRadiativeLosses();
    energy_next = m_energy - (dt_total/(double)subcycles_rad)*m_grids[rad];
    clampWallBoundaries(dummy_grid,dummy_grid,dummy_grid,energy_next);
    m_energy = energy_next.max(energy_floor);
    m_press = (GAMMA - 1.0)*(m_energy - nonthermal_energy);
    m_temp = M_I*m_press/(2.0*K_B*m_grids[rho]);
    // Enforce temp floor everywhere
    m_temp = m_temp.max(temp_min);
    // Recompute energy after flooring temperature
    m_grids[press] = 2.0*K_B*m_grids[rho]*m_temp/M_I;
    m_energy = m_press/(GAMMA - 1.0) + nonthermal_energy;
  }
}

//Computes the magnetic magnitude, direction, and pressure terms
//based on the current values of b_x, b_y, and b_z
void PlasmaDomain::computeMagneticTerms()
{
  Grid &m_b_x = m_grids[b_x], &m_b_y = m_grids[b_y], &m_b_z = m_grids[b_z];
  m_grids[b_magnitude] = (m_b_x.square() + m_b_y.square() + m_b_z.square()).sqrt();
  m_grids[b_hat_x] = m_b_x/m_grids[b_magnitude];
  m_grids[b_hat_y] = m_b_y/m_grids[b_magnitude];

  m_grids[mag_press] = (m_b_x*m_b_x + m_b_y*m_b_y + m_b_z*m_b_z)/(8.0*PI);
  m_grids[mag_pxx] = (-m_b_x*m_b_x + m_b_y*m_b_y + m_b_z*m_b_z)/(8.0*PI);
  m_grids[mag_pyy] = (m_b_x*m_b_x - m_b_y*m_b_y + m_b_z*m_b_z)/(8.0*PI);
  m_grids[mag_pzz] = (m_b_x*m_b_x + m_b_y*m_b_y - m_b_z*m_b_z)/(8.0*PI);
  m_grids[mag_pxy] = -m_b_x*m_b_y/(4.0*PI);
  m_grids[mag_pxz] = -m_b_x*m_b_z/(4.0*PI);
  m_grids[mag_pyz] = -m_b_y*m_b_z/(4.0*PI);
}

//Recompute pressure, energy, velocity, radiative loss rate, dt, dt_thermal, dt_rad
//from the current base variables. Does not compute variables that are shut off
//by the member physics settings.
void PlasmaDomain::recomputeDerivedVariables()
{
  m_grids[press] = 2.0*K_B*m_grids[rho]*m_grids[temp]/M_I;
  m_grids[energy] = m_grids[press]/(GAMMA - 1.0)
                  + 0.5*(m_grids[mom_x].square() + m_grids[mom_y].square())/m_grids[rho]
                  + m_grids[mag_press];
  m_grids[v_x] = m_grids[mom_x]/m_grids[rho];
  m_grids[v_y] = m_grids[mom_y]/m_grids[rho];
  recomputeRadiativeLosses();
  recomputeDT();
  if(thermal_conduction) recomputeDTThermal();
  if(radiative_losses) recomputeDTRadiative();
}

void PlasmaDomain::recomputeTemperature()
{
  Grid &m_press = m_grids[press];
  m_press = (GAMMA - 1.0)*(m_grids[energy] - 0.5*(m_grids[mom_x].square() + m_grids[mom_y].square())/m_grids[rho] - m_grids[mag_press]);
  m_press = m_press.max((GAMMA - 1.0)*thermal_energy_min);
  m_grids[temp] = (M_I*m_press/(2.0*K_B*m_grids[rho])).max(temp_min);
}

void PlasmaDomain::catchUnderdensity()
{
  for(int i=0; i<m_xdim; i++){
    for(int j=0; j<m_ydim; j++){
      if(m_grids[rho](i,j) < rho_min){
        std::cout << "Density too low at (" << i << "," << j << ")" << "\n";
        m_grids[rho](i,j) = rho_min;
      }
    }
  }
}

//Enforces zero motion and unchanging rho and energy at Wall boundaries
void PlasmaDomain::clampWallBoundaries(Grid& mom_x_next, Grid& mom_y_next, Grid& rho_next, Grid& energy_next)
{
  if(y_bound_1 == BoundaryCondition::Wall || y_bound_2 == BoundaryCondition::Wall){
    if(y_bound_1 == BoundaryCondition::Wall) for(int i=0; i<m_xdim; i++){
      mom_x_next(i,0) = 0.0;
      mom_y_next(i,0) = 0.0;
      rho_next(i,0) = m_grids[rho](i,0);
      energy_next(i,0) = m_grids[energy](i,0);
    }
    if(y_bound_2 == BoundaryCondition::Wall) for(int i=0; i<m_xdim; i++){
      mom_x_next(i,m_ydim-1) = 0.0;
      mom_y_next(i,m_ydim-1) = 0.0;
      rho_next(i,m_ydim-1) = m_grids[rho](i,m_ydim-1);
      energy_next(i,m_ydim-1) = m_grids[energy](i,m_ydim-1);
    }
  }
  if(x_bound_1 == BoundaryCondition::Wall || x_bound_2 == BoundaryCondition::Wall){
    if(x_bound_1 == BoundaryCondition::Wall) for(int j=0; j<m_ydim; j++){
      mom_x_next(0,j) = 0.0;
      mom_y_next(0,j) = 0.0;
      rho_next(0,j) = m_grids[rho](0,j);
      energy_next(0,j) = m_grids[energy](0,j);
    }
    if(x_bound_2 == BoundaryCondition::Wall) for(int j=0; j<m_ydim; j++){
      mom_x_next(m_xdim-1,j) = 0.0;
      mom_y_next(m_xdim-1,j) = 0.0;
      rho_next(m_xdim-1,j) = m_grids[rho](m_xdim-1,j);
      energy_next(m_xdim-1,j) = m_grids[energy](m_xdim-1,j);
    }
  }
}

void PlasmaDomain::recomputeDT()
{
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid c_s = (GAMMA*m_grids[press]/m_grids[rho]).sqrt();
  Grid abs_vx = m_dx/(c_s+m_grids[v_x].abs());
  Grid abs_vy = m_dy/(c_s+m_grids[v_y].abs());
  m_grids[dt] = abs_vx.min(abs_vy);
}

void PlasmaDomain::recomputeDTThermal()
{
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  if(!flux_saturation){
    m_grids[dt_thermal] = K_B/KAPPA_0*(m_grids[rho]/M_I)*m_dx*m_dy/m_grids[temp].pow(2.5);
  } else {
    Grid kappa_modified(m_xdim,m_ydim,0.0);
    Grid con_flux_x(m_xdim,m_ydim,0.0);
    Grid con_flux_y(m_xdim,m_ydim,0.0);
    Grid field_temp_gradient = derivative1D(m_grids[temp],0)*m_grids[b_hat_x]
                              + derivative1D(m_grids[temp],1)*m_grids[b_hat_y];
    fieldAlignedConductiveFlux(con_flux_x, con_flux_y, m_grids[temp], m_grids[rho],
                                  m_grids[b_hat_x], m_grids[b_hat_y], KAPPA_0);
    Grid flux_c = (con_flux_x.square() + con_flux_y.square()).sqrt();
    saturateConductiveFlux(con_flux_x, con_flux_y, m_grids[rho], m_grids[temp]);
    Grid flux_sat = (con_flux_x.square() + con_flux_y.square()).sqrt();
    kappa_modified = (flux_sat/field_temp_gradient).abs();
    m_grids[dt_thermal] = K_B/kappa_modified*(m_grids[rho]/M_I)*m_dx*m_dy;
  }
}

void PlasmaDomain::recomputeDTRadiative()
{
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  m_grids[dt_rad] = (m_grids[energy]/m_grids[rad]).abs();
}


//Compute the radiative loss rate (stored in m_grids[rad])
//using piecewise approximation for optically thin coronal radiative losses
void PlasmaDomain::recomputeRadiativeLosses()
{
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  #pragma omp parallel
  {
    #if BENCHMARKING_ON
    InstrumentationTimer timer("loop thread");
    #endif
    #pragma omp for collapse(2)
    for(int i=0; i<m_xdim; i++){
      for(int j=0; j<m_ydim; j++){
        if(m_grids[temp](i,j) < temp_chromosphere){
          m_grids[rad](i,j) = 0.0;
        }
        else {
          double logtemp = std::log10(m_grids[temp](i,j));
          double n = m_grids[rho](i,j)/M_I;
          double chi, alpha;
          if(logtemp <= 4.97){
            chi = 1.09e-31; //also adjust chi to ensure continuity
            alpha = 2.0; //alpha 3 might be better approx?
            // chi = 1.17e-36;
            // alpha = 3.0;
          } else if(logtemp <= 5.67){
            chi = 8.87e-17;
            alpha = -1.0;
          } else if(logtemp <= 6.18){
            chi = 1.90e-22;
            alpha = 0.0;
          } else if(logtemp <= 6.55){
            chi = 3.53e-13;
            alpha = -1.5;
          } else if(logtemp <= 6.90){
            chi = 3.46e-25;
            alpha = 1.0/3.0;
          } else if(logtemp <= 7.63){
            chi = 5.49e-16;
            alpha = -1.0;
          } else{
            chi = 1.96e-27;
            alpha = 0.5;
          }
          m_grids[rad](i,j) = n*n*chi*std::pow(m_grids[temp](i,j),alpha);
          if(m_grids[temp](i,j) < temp_chromosphere + radiation_ramp){
            double ramp = (m_grids[temp](i,j) - temp_chromosphere)/radiation_ramp;
            m_grids[rad](i,j) *= ramp;
          }
        }
      }
    }
  }
}

//Computes 1D cell-centered conductive flux from temperature "temp"
//Flux computed in direction indicated by "index": 0 for x, 1 for y
//k0 is conductive coefficient
Grid PlasmaDomain::oneDimConductiveFlux(const Grid &temp, const Grid &rho, double k0, int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid kappa_max = m_dx*m_dy*K_B*(rho/M_I)/dt_thermal_min;
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
void PlasmaDomain::fieldAlignedConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &temp, const Grid &rho,
                                    const Grid &b_hat_x, const Grid &b_hat_y, double k0){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  int xdim = temp.rows();
  int ydim = temp.cols();
  Grid con_flux_x = oneDimConductiveFlux(temp, rho, k0, 0);
  Grid con_flux_y = oneDimConductiveFlux(temp, rho, k0, 1);
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
void PlasmaDomain::saturateConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &rho, const Grid &temp){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid sat_flux_mag = (1.0/6.0)*(3.0/2.0)*(rho/M_I)*(K_B*temp).pow(1.5)/std::sqrt(M_ELECTRON);
  Grid flux_mag = (flux_out_x.square() + flux_out_y.square()).sqrt();
  Grid scale_factor = sat_flux_mag /((sat_flux_mag.square() + flux_mag.square()).sqrt());
  flux_out_x *= scale_factor;
  flux_out_y *= scale_factor;
}

