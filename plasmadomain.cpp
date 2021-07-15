#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include "utils.hpp"
#include "derivs.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"

//Enum allows each variable's index to be referred to directly
enum Variable{
  rho,temp,mom_x,mom_y,b_x,b_y,b_z, //base variables (carried over between iterations)
  press,energy,rad,dt,dt_thermal,dt_rad,v_x,v_y, //derived variables (derived from base variables)
  grav_x,grav_y,b_magnitude,b_hat_x,b_hat_y, //constant variables (unchanging bw iterations)
  mag_press,mag_pxx,mag_pyy,mag_pzz,mag_pxy,mag_pxz,mag_pyz, //constant variables
  num_variables //never add new variable after this in the enum!
};

//Corresponding variable names for file I/O
//Must match the ordering of the Variable enum above
const std::vector<std::string> PlasmaDomain::m_var_names = {
  "rho","temp","mom_x","mom_y","b_x","b_y","b_z",
  "press","energy","rad","dt","dt_thermal","dt_rad","v_x","v_y",
  "grav_x","grav_y","b_magnitude","b_hat_x","b_hat_y",
  "mag_press","mag_pxx","mag_pyy","mag_pzz","mag_pxy","mag_pxz","mag_pyz"
};

PlasmaDomain::PlasmaDomain() : PlasmaDomain::PlasmaDomain(1,1,1.0,1.0,"run") {}
PlasmaDomain::PlasmaDomain(const char* run_name) : PlasmaDomain::PlasmaDomain(1,1,1.0,1.0,run_name) {}

PlasmaDomain::PlasmaDomain(size_t xdim, size_t ydim, double dx, double dy, const char* run_name)
{
  m_xdim = xdim; m_ydim = ydim;
  m_dx = dx; m_dy = dy;
  m_time = 0.0; m_iter = 0;
  m_run_name = std::string(run_name);
  m_out_file.open(m_run_name+".out");
  for(int i=0; i<num_variables; i++){
    m_grids.push_back(Grid(xdim,ydim,0.0));
    m_output_flags.push_back(false);
    m_state_flags.push_back(false);
  }
}

void PlasmaDomain::hydrostaticInitialize()
{
  double base_rho = M_I*1.0e12; //initial mass density at base, g cm^-3
  double scale_height = 2.0*K_B*temp_chromosphere/(M_I*BASE_GRAV);
  m_grids[rho] = HydrostaticFalloff(base_rho,scale_height,m_xdim,m_ydim);
  m_grids[temp] = Grid(m_xdim,m_ydim,temp_chromosphere);
  m_grids[b_x] = BipolarField(m_xdim, m_ydim, b_0, scale_height, 0);
  m_grids[b_y] = BipolarField(m_xdim, m_ydim, b_0, scale_height, 1);
  computeMagneticTerms();
  recomputeDerivedVariables();
}

void PlasmaDomain::gaussianInitialize()
{
  //Simple Gaussian test case (no magnetic field)
  m_grids[rho] = GaussianGrid(m_xdim, m_ydim, 1.0e-15, 1.0e-10);
  m_grids[mom_x] = Grid::Zero(m_xdim,m_ydim);
  m_grids[mom_y] = Grid::Zero(m_xdim,m_ydim);
  m_grids[temp] = GaussianGrid(m_xdim, m_ydim, temp_chromosphere, 2.0*temp_chromosphere);
  computeMagneticTerms();
  recomputeDerivedVariables();
  // double b0 = 0.1;
  // double h = (m_xdim*DX)/(4.0*PI);
  // b_x = BipolarField(m_xdim, m_ydim, b0, h, 0);
  // b_y = BipolarField(m_xdim, m_ydim, b0, h, 1);
  // b_z = Grid::Zero(m_xdim,m_ydim);
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

  if(verbose){
    printf("Iter: %i/%i|t: %f|dt: %f",m_iter,n_iterations,m_time,min_dt);
    if(thermal_conduction) printf("|Cond. Subcyc: %i",subcycles_thermal);
    if(radiative_losses) printf("|Rad. Subcyc: %i",subcycles_rad);
    printf("\n");
  }
  if(thermal_conduction) subcycleConduction(subcycles_thermal,min_dt);
  if(radiative_losses) subcycleRadiation(subcycles_rad,min_dt);

  //Advance time by min_dt
  Grid &m_mom_x = m_grids[mom_x], &m_mom_y = m_grids[mom_y],
        &m_v_x = m_grids[v_x], &m_v_y = m_grids[v_y],
        &m_rho = m_grids[rho], &m_energy = m_grids[energy], &m_press = m_grids[press];
  Grid viscous_force_x = epsilon_viscous*(0.5*m_dx*m_dx/dt_raw)*laplacian(m_mom_x);
  Grid viscous_force_y = epsilon_viscous*(0.5*m_dy*m_dy/dt_raw)*laplacian(m_mom_y);
  Grid zero(1,1,0.0);
  Grid rho_next = m_rho - min_dt*divergence(m_rho,zero,zero,m_v_x,m_v_y);
  Grid mom_x_next = m_mom_x - min_dt*divergence(m_mom_x, m_press + m_grids[mag_pxx], m_grids[mag_pxy], m_v_x, m_v_y)
                    + min_dt*m_rho*m_grids[grav_x] + min_dt*viscous_force_x;
  Grid mom_y_next = m_mom_y - min_dt*divergence(m_mom_y, m_grids[mag_pxy], m_press + m_grids[mag_pyy], m_v_x, m_v_y)
                    + min_dt*m_rho*m_grids[grav_y] + min_dt*viscous_force_y;
  Grid energy_next = m_energy - min_dt*divergence(m_energy+m_press, m_grids[mag_pxx]*m_v_x + m_grids[mag_pxy]*m_v_y, m_grids[mag_pxy]*m_v_x + m_grids[mag_pyy]*m_v_y, m_v_x, m_v_y) 
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
  for(int subcycle = 0; subcycle < subcycles_thermal; subcycle++){
    Grid con_flux_x(m_xdim,m_ydim,0.0);
    Grid con_flux_y(m_xdim,m_ydim,0.0);
    field_aligned_conductive_flux(con_flux_x, con_flux_y, m_temp, m_grids[rho], m_grids[b_hat_x], m_grids[b_hat_y], KAPPA_0);
    if(flux_saturation) saturate_conductive_flux(con_flux_x, con_flux_y, m_grids[rho], m_temp);
    m_energy = m_energy - (dt_total/(double)subcycles_thermal)*(derivative1D(con_flux_x,0)+derivative1D(con_flux_y,1));
    m_energy = m_energy.max(energy_floor);
    m_press = (GAMMA - 1.0)*(m_energy - nonthermal_energy);
    m_temp = M_I*m_press/(2.0*K_B*m_grids[rho]);
    // Enforce constant chromospheric temperature
    for(int i=0; i<m_xdim; i++) for(int j=0; j<chromosphere_depth; j++) m_temp(i,j) = temp_chromosphere;
    // Enforce temp floor everywhere
    m_temp = m_temp.max(temp_chromosphere);
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
  for(int subcycle = 0; subcycle < subcycles_rad; subcycle++){
    recomputeRadiativeLosses();
    m_energy = m_energy - (dt_total/(double)subcycles_rad)*m_grids[rad];
    m_energy = m_energy.max(energy_floor);
    m_press = (GAMMA - 1.0)*(m_energy - nonthermal_energy);
    m_temp = M_I*m_press/(2.0*K_B*m_grids[rho]);
    // Enforce constant chromospheric temperature
    // for(int i=0; i<m_xdim; i++) for(int j=0; j<chromosphere_depth; j++) m_temp(i,j) = temp_chromosphere;
    // Enforce temp floor everywhere
    m_temp = m_temp.max(temp_chromosphere);
    // Recompute energy after flooring temperature
    m_grids[press] = 2.0*K_B*m_grids[rho]*m_temp/M_I;
    m_energy = m_press/(GAMMA - 1.0) + nonthermal_energy;
  }
}

//Read in variables from .state file
//This function will abort execution if an invalid variable name is encountered
//Does not check that the variables read from the file are a complete
//description of the plasma, nor that none of them are contradictory
void PlasmaDomain::readStateFile(const char* in_filename)
{
  std::ifstream in_file(in_filename);

  //Read in grid dimensions
  std::string line, element;
  std::getline(in_file, line);
  std::istringstream ss_dim(line);
  std::getline(ss_dim,element,',');
  m_xdim = std::stoi(element);
  std::getline(ss_dim,element,',');
  m_ydim = std::stoi(element);

  //Read in physical dimensions
  std::getline(in_file,line);
  std::istringstream ss_dxdy(line);
  std::getline(ss_dxdy,element,',');
  m_dx = std::stod(element);
  std::getline(ss_dxdy,element,',');
  m_dy = std::stod(element);

  //Read in time of state file
  std::getline(in_file,line);
  std::istringstream ss_time(line);
  std::getline(ss_time,element,'=');
  assert(element == "t");
  std::getline(ss_time,element);
  m_time = std::stod(element);

  while(std::getline(in_file, line)){
    //Ensure that variable in file is valid
    auto it = std::find(m_var_names.begin(),m_var_names.end(),line);
    assert(it != m_var_names.end());
    auto index = std::distance(m_var_names.begin(),it);
    Grid curr_grid(m_xdim,m_ydim);
    //Read in Grid corresponding to variable
    int i=0; int j=0;
    std::string row; std::string el;
    std::getline(in_file,line);
    i=0;
    std::istringstream ss_line(line);
    while(std::getline(ss_line,row,';')){
      j=0;
      std::istringstream ss_row(row);
      while(std::getline(ss_row,el,',')){
        curr_grid(i,j) = std::stod(el);
        j++;
      }
      i++;
    }
    m_grids[index] = curr_grid;
  }
  in_file.close();
  computeMagneticTerms();
  recomputeDerivedVariables();
}

void PlasmaDomain::setDefaultSettings()
{
  x_bound_1 = static_cast<BoundaryCondition>(XBOUND1);
  x_bound_2 = static_cast<BoundaryCondition>(XBOUND2);
  y_bound_1 = static_cast<BoundaryCondition>(YBOUND1);
  y_bound_2 = static_cast<BoundaryCondition>(YBOUND2);
  radiative_losses = RADIATIVE_LOSSES_ON;
  ambient_heating = AMBIENT_HEATING_ON;
  thermal_conduction = THERMAL_CONDUCTION_ON;
  flux_saturation = FLUX_SATURATION;
  temp_chromosphere = TEMP_CHROMOSPHERE; //Minimum allowed temperature
  radiation_ramp = RADIATION_RAMP;  //Width of cutoff ramp, in units of temperature, for low-temp radiation
  heating_rate = HEATING_RATE;  //Constant ambient heating rate
  b_0 = B0;  //Base value of magnetic field
  chromosphere_depth = CHROMOSPHERE_DEPTH; //In number of Grid cells
  epsilon = EPSILON; epsilon_thermal = EPSILON_THERMAL; epsilon_rad = EPSILON_RADIATIVE; //Time step calculation
  epsilon_viscous = EPSILON_VISCOUS; //Prefactor for artificial viscosity
  dt_thermal_min = DT_THERMAL_MIN; //Minimum timestep for thermal conduction
  rho_min = GRIDFLOOR;
  thermal_energy_min = THERMALENERGYFLOOR; //Lower bounds for mass density and thermal energy density
  n_iterations = NT;
  output_interval = OUTPUT_INTERVAL;
  m_output_flags[rho] = RHO_OUT;
  m_output_flags[temp] = TEMP_OUT;
  m_output_flags[press] = PRESS_OUT;
  m_output_flags[rad] = RAD_OUT;
  m_output_flags[energy] = ENERGY_OUT;
  m_output_flags[v_x] = VX_OUT;
  m_output_flags[v_y] = VY_OUT;
  m_output_flags[dt] = DT_OUT;
  m_output_flags[dt_thermal] = DT_THERMAL_OUT;
  m_output_flags[dt_rad] = DT_RAD_OUT;
  m_state_flags[rho] = true;
  m_state_flags[mom_x] = true;
  m_state_flags[mom_y] = true;
  m_state_flags[temp] = true;
  m_state_flags[b_x] = true;
  m_state_flags[b_y] = true;
  m_state_flags[b_z] = true;
}

void PlasmaDomain::outputPreamble()
{
  m_out_file << m_xdim << "," << m_ydim << std::endl;
  m_out_file << m_dx << "," << m_dy << std::endl;
  m_out_file << m_grids[b_x].format(',',';');
  m_out_file << m_grids[b_y].format(',',';');
}

void PlasmaDomain::outputCurrentState()
{
  m_out_file << "t=" << m_time << std::endl;
  for(int i=0; i<num_variables; i++){
    if(m_output_flags[i]){
      m_out_file << m_var_names[i] << std::endl;
      m_out_file << m_grids[i].format(',',';');
    }
  }
}

//Output current state into state file (toggles between overwriting two different files
//so that the most recent state is still written in event of a crash)
//All variables marked in m_state_flags are output here
void PlasmaDomain::writeStateFile() const
{
  std::ofstream state_file;
  state_file.open(m_run_name+std::to_string(m_iter%2)+".state");
  state_file << m_xdim << "," << m_ydim << std::endl;
  state_file << m_dx << "," << m_dy << std::endl;
  state_file << "t=" << m_time << std::endl;
  for(int i=0; i<num_variables; i++){
    if(m_state_flags[i]){
      state_file << m_var_names[i] << std::endl;
      state_file << m_grids[i].format(',',';',-1);
    }
  }
  state_file.close();
}

//Gets rid of the older of the two state files at the end of the sim run
void PlasmaDomain::cleanUpStateFiles() const
{
  rename((m_run_name+std::to_string((m_iter-1)%2)+".state").c_str(),
          (m_run_name+".state").c_str());
  remove((m_run_name+std::to_string((m_iter-2)%2)+".state").c_str());
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
  m_grids[temp] = (M_I*m_press/(2.0*K_B*m_grids[rho])).max(temp_chromosphere);
  //Enforce chromosphere temperature
  for(int i=0; i<m_xdim; i++){
    for(int j=0; j<chromosphere_depth; j++){
      m_grids[temp](i,j) = temp_chromosphere;
    }
  }
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
    field_aligned_conductive_flux(con_flux_x, con_flux_y, m_grids[temp], m_grids[rho],
                                  m_grids[b_hat_x], m_grids[b_hat_y], KAPPA_0);
    Grid flux_c = (con_flux_x.square() + con_flux_y.square()).sqrt();
    saturate_conductive_flux(con_flux_x, con_flux_y, m_grids[rho], m_grids[temp]);
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

//Sets gravity to fall off from base_gravity at bottom of the domain,
//as though from the surface of a planet with radius r_solar
void PlasmaDomain::setSolarGravity(double base_gravity, double r_solar)
{
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid& m_grav_y = m_grids[grav_y];
  for(int j=0; j<m_ydim; j++){
    double y = j*m_dy;
    for(int i=0; i<m_xdim; i++){
      m_grav_y(i,j) = base_gravity*std::pow(r_solar/(r_solar+y),2.0);
    }
  }
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
            // chi = 1.09e-31; //also adjust chi to ensure continuity
            // alpha = 2.0; //alpha 3 might be better approx?
            chi = 1.17e-36;
            alpha = 3.0;
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
            //Also try using linear ramp
            // double ramp = (m_grids[temp](i,j) - temp_chromosphere)/radiation_ramp;
            double ramp = 0.5*(1.0 - std::cos((m_grids[temp](i,j) - temp_chromosphere)*PI/radiation_ramp));
            m_grids[rad](i,j) *= ramp;
          }
        }
      }
    }
  }
}

void PlasmaDomain::setOutputFlag(int var, bool new_flag) { m_output_flags[var] = new_flag; }

void PlasmaDomain::setOutputFlags(const std::vector<int> vars, bool new_flag)
{
  for(int var : vars){
    m_output_flags[var] = new_flag;
  }
}

void PlasmaDomain::setStateFlag(int var, bool new_flag) { m_state_flags[var] = new_flag; }

void PlasmaDomain::setStateFlags(const std::vector<int> vars, bool new_flag)
{
  for(int var : vars){
    m_state_flags[var] = new_flag;
  }
}
