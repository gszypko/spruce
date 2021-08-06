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

PlasmaDomain::PlasmaDomain(const char* run_name, const char* config_filename)
{
  m_xdim = 1; m_ydim = 1; //default values, overwritten by initialize()
  m_run_name = std::string(run_name);
  for(int i=0; i<num_variables; i++){
    m_output_flags.push_back(false);
    m_state_flags.push_back(false);
  }
  readConfigFile(config_filename);
  m_out_file.open(m_run_name+".out");
  m_time = 0.0; m_iter = 0;
  for(int i=0; i<num_variables; i++) m_grids.push_back(Grid(m_xdim,m_ydim,0.0));
  setStateFlags({"rho","mom_x","mom_y","temp","b_x","b_y","b_z","pos_x","pos_y","grav_x","grav_y"},true);
}

void PlasmaDomain::initialize(const Grid& rho, const Grid& temp, const Grid& mom_x, const Grid& mom_y,
                              const Grid& b_x, const Grid& b_y, const Grid& b_z,
                              const Grid& pos_x, const Grid& pos_y, const Grid& grav_x, const Grid& grav_y)
{
  m_xdim = rho.rows(); m_ydim = rho.cols();
  for(int var=0; var<num_variables; var++){
    m_grids[var] = Grid(m_xdim,m_ydim,0.0);
  }
  m_grids[Variable::rho] = rho;
  m_grids[Variable::temp] = temp;
  m_grids[Variable::mom_x] = mom_x;
  m_grids[Variable::mom_y] = mom_y;
  m_grids[Variable::b_x] = b_x;
  m_grids[Variable::b_y] = b_y;
  m_grids[Variable::b_z] = b_z;
  m_grids[Variable::pos_x] = pos_x;
  m_grids[Variable::pos_y] = pos_y;
  m_grids[Variable::grav_x] = grav_x;
  m_grids[Variable::grav_y] = grav_y;
  computeIterationBounds();
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
}

//Compute lower and upper x- and y- indicies for differential operations
//from boundary conditions, to exclude ghost zones. Results stored in m_xl, m_xu, m_yl, m_yu
void PlasmaDomain::computeIterationBounds()
{
  assert(m_xdim > 2*N_GHOST && m_ydim > 2*N_GHOST && "Grid too small for ghost zones");
  if(x_bound_1 == BoundaryCondition::Open || x_bound_1 == BoundaryCondition::Wall) m_xl = N_GHOST;
  else{ assert(x_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xl = 0; }
  if(x_bound_2 == BoundaryCondition::Open || x_bound_2 == BoundaryCondition::Wall) m_xu = m_xdim - N_GHOST - 1;
  else{ assert(x_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xu = m_xdim - 1; }
  if(y_bound_1 == BoundaryCondition::Open || y_bound_1 == BoundaryCondition::Wall) m_yl = N_GHOST;
  else{ assert(y_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yl = 0; }
  if(y_bound_2 == BoundaryCondition::Open || y_bound_2 == BoundaryCondition::Wall) m_yu = m_ydim - N_GHOST - 1;
  else{ assert(y_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yu = m_ydim - 1; }
}