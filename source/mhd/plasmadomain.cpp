//plasmadomain.cpp
//PlasmaDomain functionality related to initializing
//and setting/modifying simulation configuration

#include <vector>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include "utils.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"

PlasmaDomain::PlasmaDomain() : PlasmaDomain::PlasmaDomain(1,1,1.0,1.0,"run") {}
PlasmaDomain::PlasmaDomain(const char* run_name) : PlasmaDomain::PlasmaDomain(1,1,1.0,1.0,run_name) {}

PlasmaDomain::PlasmaDomain(const char* run_name, const char* config_filename)
{
  m_run_name = std::string(run_name);
  for(int i=0; i<num_variables; i++){
    m_output_flags.push_back(false);
    m_state_flags.push_back(false);
  }
  readConfigFile(config_filename);
  computeIterationBounds();
  m_out_file.open(m_run_name+".out");
  m_time = 0.0; m_iter = 0;
  for(int i=0; i<num_variables; i++) m_grids.push_back(Grid(m_xdim,m_ydim,0.0));
}

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
  computeIterationBounds();
}

void PlasmaDomain::hydrostaticInitialize()
{
  double base_rho = M_I*1.0e12; //initial mass density at base, g cm^-3
  double scale_height = 2.0*K_B*temp_chromosphere/(M_I*BASE_GRAV);
  Grid &m_pos_x = m_grids[pos_x], &m_pos_y = m_grids[pos_y];
  for(int i=0; i<m_xdim; i++){
    for(int j=0; j<m_ydim; j++){
      m_pos_x(i,j) = (i - (double)(m_xdim-1)*0.5)*m_dx;
      m_pos_y(i,j) = j*m_dy;
    }
  }
  m_grids[rho] = HydrostaticFalloff(base_rho,scale_height,m_pos_y);
  m_grids[temp] = Grid(m_xdim,m_ydim,temp_chromosphere);
  m_grids[b_x] = BipolarField(m_pos_x, m_pos_y, b_0, scale_height, 0);
  m_grids[b_y] = BipolarField(m_pos_x, m_pos_y, b_0, scale_height, 1);
  computeConstantTerms();
  recomputeDerivedVariables();
}

//Generates gaussian distribution in density and temperature, with given min
//and max density and temperature, distributed with the given standard deviations
//in x and y (in units of grid cell widths)
void PlasmaDomain::gaussianInitialize(double min_rho, double max_rho, double min_temp, double max_temp,
                                      double std_dev_x, double std_dev_y)
{
  m_grids[rho] = GaussianGrid(m_xdim, m_ydim, min_rho, max_rho, std_dev_x, std_dev_y);
  m_grids[mom_x] = Grid::Zero(m_xdim,m_ydim);
  m_grids[mom_y] = Grid::Zero(m_xdim,m_ydim);
  m_grids[temp] = GaussianGrid(m_xdim, m_ydim, min_temp, max_temp, std_dev_x, std_dev_y);
  computeConstantTerms();
  recomputeDerivedVariables();
  // double b0 = 0.1;
  // double h = (m_xdim*DX)/(4.0*PI);
  // b_x = BipolarField(m_xdim, m_ydim, b0, h, 0);
  // b_y = BipolarField(m_xdim, m_ydim, b0, h, 1);
  // b_z = Grid::Zero(m_xdim,m_ydim);
}

//Sets gravity to fall off from base_gravity at bottom of the domain,
//as though from the surface of a planet with radius r_solar
void PlasmaDomain::setSolarGravity(double base_gravity, double r_solar)
{
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid &m_grav_y = m_grids[grav_y], &m_pos_y = m_grids[pos_y];
  for(int j=0; j<m_ydim; j++){
    for(int i=0; i<m_xdim; i++){
      m_grav_y(i,j) = -base_gravity*std::pow(r_solar/(r_solar+m_pos_y(i,j)),2.0);
    }
  }
}

//Compute lower and upper x- and y- indicies for differential operations
//from boundary conditions, to exclude ghost zones. Results stored in m_xl, m_xu, m_yl, m_yu
void PlasmaDomain::computeIterationBounds()
{
  if(x_bound_1 == BoundaryCondition::Open || x_bound_1 == BoundaryCondition::Wall) m_xl = N_GHOST;
  else{ assert(x_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xl = 0; }
  if(x_bound_2 == BoundaryCondition::Open || x_bound_2 == BoundaryCondition::Wall) m_xu = m_xdim - N_GHOST - 1;
  else{ assert(x_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xu = m_xdim - 1; }
  if(y_bound_1 == BoundaryCondition::Open || y_bound_1 == BoundaryCondition::Wall) m_yl = N_GHOST;
  else{ assert(y_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yl = 0; }
  if(y_bound_2 == BoundaryCondition::Open || y_bound_2 == BoundaryCondition::Wall) m_yu = m_ydim - N_GHOST - 1;
  else{ assert(y_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yu = m_ydim - 1; }
}