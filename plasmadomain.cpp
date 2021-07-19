//plasmadomain.cpp
//PlasmaDomain functionality related to initializing
//and setting/modifying simulation settings

#include <vector>
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
  epsilon = EPSILON; epsilon_thermal = EPSILON_THERMAL; epsilon_rad = EPSILON_RADIATIVE; //Time step calculation
  epsilon_viscous = EPSILON_VISCOUS; //Prefactor for artificial viscosity
  dt_thermal_min = DT_THERMAL_MIN; //Minimum timestep for thermal conduction
  rho_min = RHO_MIN;
  temp_min = TEMP_MIN;
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
