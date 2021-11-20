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
namespace fs = std::filesystem;
#include "utils.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"

PlasmaDomain::PlasmaDomain(const fs::path &out_path, const fs::path &config_path, bool continue_mode)
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
  m_out_directory = out_path;
  fs::path out_filename("mhd.out");
  if(!continue_mode){
    fs::copy(config_path, m_out_directory/(config_path.filename()), fs::copy_options::overwrite_existing);
    m_out_file.open(m_out_directory/out_filename);
  } else {
    m_out_file.open(m_out_directory/out_filename,std::ofstream::app);
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
  for(int i : state_vars)
    m_grids[i] = input_vars[i];
  computeIterationBounds();
  computeConstantTerms();
  recomputeDerivedVariables();
  writeStateFile("init");
}

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
  if(origin_position.compare("center") == 0){
    pos -= 0.5*(pos(xdim-1,ydim-1) + 0.5*d(xdim-1,ydim-1));
  } else if(origin_position.compare("upper") == 0){
    pos -= (pos(xdim-1,ydim-1) + 0.5*d(xdim-1,ydim-1));
  } else assert(origin_position == "lower" && "origin_position must be one of upper, lower, or center");

  return pos;
}