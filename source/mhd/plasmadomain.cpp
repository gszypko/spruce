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

//Construction with newly constructed initial state (not continue mode)
PlasmaDomain::PlasmaDomain(const fs::path &out_path, const fs::path &config_path, const std::vector<Grid>& input_vars,
                          double ion_mass, double adiabatic_index) : m_module_handler(*this)
{
  m_ion_mass = ion_mass;
  m_adiabatic_index = adiabatic_index;
  m_xdim = input_vars[0].rows(); m_ydim = input_vars[0].cols();
  for(int var=0; var<num_variables; var++){
    m_grids.push_back(Grid(m_xdim,m_ydim,0.0));
  }
  for(int i : state_vars)
    m_grids[i] = input_vars[i];
  m_time = 0.0; 
  configureSimulation(config_path,out_path,false);
}

//Construction with initial state from state file (continue mode)
PlasmaDomain::PlasmaDomain(const fs::path &out_path, const fs::path &config_path, const fs::path &state_file, bool continue_mode) : m_module_handler(*this)
{
  readStateFile(state_file,continue_mode);
  configureSimulation(config_path,out_path,continue_mode);
}

void PlasmaDomain::configureSimulation(const fs::path &config_path, const fs::path &out_path, bool continue_mode)
{
  //DEFAULT VALUES, TO BE OVERWRITTEN BY readConfigFile
  open_boundary_strength = 1.0;
  time_integrator = TimeIntegrator::Euler;
  std_out_interval = 1;
  safe_state_mode = true;
  //***************************************************
  m_iter = 0;
  this->continue_mode = continue_mode;
  for(int i=0; i<num_variables; i++){
    m_output_flags.push_back(false);
  }
  m_out_directory = out_path;
  fs::path out_filename("mhd.out");
  if(continue_mode){
    m_out_file.open(m_out_directory/out_filename,std::ofstream::app);
  } else {
    m_out_file.open(m_out_directory/out_filename);
    fs::path new_config_path = m_out_directory/(config_path.filename());
    fs::directory_entry new_config_dir(new_config_path);
    if(!(new_config_dir.exists() && fs::equivalent(config_path,new_config_path))){
      std::cout << "Copying " << config_path.string() << " into " << m_out_directory.string() << "...\n";
      fs::copy(config_path, new_config_path, fs::copy_options::overwrite_existing);
    }
    else std::cout << config_path.string() << " already located in output directory.\n";
  }
  assert(validateCellSizesAndPositions(m_grids[d_x],m_grids[pos_x],0) && validateCellSizesAndPositions(m_grids[d_y],m_grids[pos_y],1) && "Cell sizes and positions must correspond");
  readConfigFile(config_path);
  computeIterationBounds();
  computeConstantTerms();
  recomputeDerivedVariables();
  if(!continue_mode) writeStateFile("init");
}

//Compute lower and upper x- and y- indicies for differential operations
//from boundary conditions, to exclude ghost zones. Results stored in m_xl, m_xu, m_yl, m_yu
//Also generates m_ghost_zone_mask Grid.
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
  m_ghost_zone_mask = Grid(m_xdim,m_ydim,0.0);
  for(int i = m_xl; i <= m_xu; i++){
    for(int j = m_yl; j <= m_yu; j++){
      m_ghost_zone_mask(i,j) = 1.0;
    }
  }
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

//Verifies that spacings d (in x-direction for index==0 and y-direction for index==1) correctly
//matches positions pos (corresponding to index as above), within fractional tolerance
bool PlasmaDomain::validateCellSizesAndPositions(const Grid& d, const Grid& pos, int index, double tolerance)
{
  int xdim = d.rows(), ydim = d.cols();
  bool valid = true;

  for(int i=0; i<xdim-1; i++){
    for(int j=0; j<ydim-1; j++){
      int i_next = i+(1-index);
      int j_next = j+index;
      double expected = pos(i,j) + 0.5*d(i,j) + 0.5*d(i_next,j_next);
      valid = valid && (std::abs((expected - pos(i_next,j_next))/expected) < tolerance);
    }
  }
  return valid;
}
