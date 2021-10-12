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
#include <filesystem>
#include "utils.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"

PlasmaDomain::PlasmaDomain(const char* out_path, const char* config_path, bool continue_mode)
{
  //DEFAULT VALUES, TO BE OVERWRITTEN BY readConfigFile
  m_xdim = 1; m_ydim = 1; open_boundary_strength = 1.0;
  std_out_interval = 1;
  safe_state_mode = true;
  //***************************************************
  m_time = 0.0; m_iter = 0;
  this->continue_mode = continue_mode;
  for(int i=0; i<num_variables; i++){
    m_output_flags.push_back(false);
  }
  readConfigFile(config_path);
  m_out_directory = std::filesystem::path(out_path);
  std::filesystem::path config_filename = std::filesystem::path(config_path).filename();
  if(!continue_mode){
    std::filesystem::copy(config_path, m_out_directory/config_filename, std::filesystem::copy_options::overwrite_existing);
    m_out_file.open(m_out_directory/"mhd.out");
  } else {
    m_out_file.open(m_out_directory/"mhd.out",std::ofstream::app);
  }
  for(int i=0; i<num_variables; i++) m_grids.push_back(Grid(m_xdim,m_ydim,0.0));
}

void PlasmaDomain::initialize(const std::vector<Grid>& input_vars, double ion_mass, double adiabatic_index)
{
  m_ion_mass = ion_mass;
  m_adiabatic_index = adiabatic_index;
  m_xdim = input_vars[0].rows(); m_ydim = input_vars[0].cols();
  for(int var=0; var<num_variables; var++){
    m_grids[var] = Grid(m_xdim,m_ydim,0.0);
  }
  for(int i=state_var_start; i < state_var_end; i++)
    m_grids[i] = input_vars[i];
  // m_grids[d_x] = input_vars[d_x];
  // m_grids[d_y] = input_vars[d_y];
  // m_grids[rho] = input_vars[rho];
  // m_grids[temp] = input_vars[temp];
  // m_grids[mom_x] = input_vars[mom_x];
  // m_grids[mom_y] = input_vars[mom_y];
  // m_grids[b_x] = input_vars[b_x];
  // m_grids[b_y] = input_vars[b_y];
  // m_grids[b_z] = input_vars[b_z];
  // m_grids[grav_x] = input_vars[grav_x];
  // m_grids[grav_y] = input_vars[grav_y];
  computeIterationBounds();
  computeConstantTerms();
  recomputeDerivedVariables();
  writeStateFile("init");
}

// //Generates gaussian distribution in density and temperature, with given min
// //and max density and temperature, distributed with the given standard deviations
// //in x and y (in units of grid cell widths)
// void PlasmaDomain::gaussianInitialize(double min_rho, double max_rho, double min_temp, double max_temp,
//                                       double std_dev_x, double std_dev_y)
// {
//   m_grids[rho] = GaussianGrid(m_xdim, m_ydim, min_rho, max_rho, std_dev_x, std_dev_y);
//   m_grids[mom_x] = Grid::Zero(m_xdim,m_ydim);
//   m_grids[mom_y] = Grid::Zero(m_xdim,m_ydim);
//   m_grids[temp] = GaussianGrid(m_xdim, m_ydim, min_temp, max_temp, std_dev_x, std_dev_y);
//   computeConstantTerms();
//   recomputeDerivedVariables();
// }

//Compute lower and upper x- and y- indicies for differential operations
//from boundary conditions, to exclude ghost zones. Results stored in m_xl, m_xu, m_yl, m_yu
void PlasmaDomain::computeIterationBounds()
{
  assert(m_xdim > 2*N_GHOST && m_ydim > 2*N_GHOST && "Grid too small for ghost zones");
  if(x_bound_1 == BoundaryCondition::Open || x_bound_1 == BoundaryCondition::Reflect || x_bound_1 == BoundaryCondition::Fixed) m_xl = N_GHOST;
  else{ assert(x_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xl = 0; }
  if(x_bound_2 == BoundaryCondition::Open || x_bound_2 == BoundaryCondition::Reflect || x_bound_2 == BoundaryCondition::Fixed) m_xu = m_xdim - N_GHOST - 1;
  else{ assert(x_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xu = m_xdim - 1; }
  if(y_bound_1 == BoundaryCondition::Open || y_bound_1 == BoundaryCondition::Reflect || y_bound_1 == BoundaryCondition::Fixed) m_yl = N_GHOST;
  else{ assert(y_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yl = 0; }
  if(y_bound_2 == BoundaryCondition::Open || y_bound_2 == BoundaryCondition::Reflect || y_bound_2 == BoundaryCondition::Fixed) m_yu = m_ydim - N_GHOST - 1;
  else{ assert(y_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yu = m_ydim - 1; }
}

//Takes not-necessarily-uniform cell sizes d and converts to cell center positions,
//returned as a Grid
//origin_position specifies location of the origin i.e. position=0 in the domain:
//"lower" is the default, "center" and "upper" are also options for each

Grid PlasmaDomain::convertCellSizesToCellPositions(const Grid& d, int index, std::string origin_position)
{
  int xdim = d.rows(), ydim = d.cols();
  Grid pos(xdim,ydim);

  // Initial computation; first assume origin at lower left
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      int i_prev = i-(1-index);
      int j_prev = j-index;
      if((i==0 && index==0) || (j==0 && index==1)) pos(i,j) = 0.5*d(i,j);
      else pos(i,j) = pos(i_prev,j_prev) + 0.5*d(i_prev,j_prev) + 0.5*d(i,j);
    }
  }

  // Translate origin to desired position
  if(origin_position == "center"){
    pos -= 0.5*(pos(xdim-1,ydim-1) + 0.5*d(xdim-1,ydim-1));
  } else if(origin_position == "upper"){
    pos -= (pos(xdim-1,ydim-1) + 0.5*d(xdim-1,ydim-1));
  } else assert(origin_position == "lower" && "origin_position must be one of upper, lower, or center");

  return pos;
}