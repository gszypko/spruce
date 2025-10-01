//plasmadomain.cpp
//PlasmaDomain functionality related to initializing
//and setting/modifying simulation configuration

#include "plasmadomain.hpp"
#include "constants.hpp"
#include <fstream>
#include <cassert>
#include <iostream>
#include <algorithm>

//Construction with initial state from state file (continue mode and custom input)
PlasmaDomain::PlasmaDomain(const fs::path &out_path, const fs::path &config_path, const fs::path &state_file, bool continue_mode, bool overwrite_init) : 
  m_module_handler(*this), m_out_directory{out_path}, m_continue_mode{continue_mode}, m_overwrite_init{overwrite_init}
{
  // initialize grids
  m_grids = std::vector<Grid>(m_gridnames.size(),Grid::Zero(1,1));

  // handle config path - if config file not within output directory, copy config file into output directory
  fs::path new_config_path = m_out_directory/(config_path.filename());
  fs::directory_entry new_config_dir(new_config_path);
  if(!fs::equivalent(config_path,new_config_path)){
    #if VERBOSE
    std::cout << "Copying " << config_path.string() << " into " << m_out_directory.string() << "...\n";
    #endif
    if(new_config_dir.exists()) fs::remove(new_config_path);
    fs::copy(config_path, new_config_path, fs::copy_options::overwrite_existing);
  }
  else{
    #if VERBOSE
    std::cout << config_path.string() << " already located in output directory.\n";
    #endif
  }
  // read from config and state files - config must be read first so equation set can be loaded
  #if VERBOSE
  std::cout << "Reading config file...\n";
  #endif
  readConfigFile(config_path);
  #if VERBOSE
  std::cout << "Reading state file...\n";
  #endif
  readStateFile(state_file,continue_mode);

  // process information loaded from config and state files
  #if VERBOSE
  std::cout << "Validating input data...\n";
  #endif
  assert(allInternalGridsInitialized() && "All internal grid quantities for PlasmaDomain must be initialized");
  assert(validateCellSizesAndPositions(m_grids[d_x],m_grids[pos_x],0) && validateCellSizesAndPositions(m_grids[d_y],m_grids[pos_y],1) && "Cell sizes and positions must correspond");
  computeIterationBounds();
  m_sg = SavitzkyGolay(m_sg_opt,Grid::Zero(m_xdim,m_ydim));
  m_eqs->setupEquationSet();
  m_module_handler.setupModules();

  if(m_multispecies_mode) {
    m_cumulative_electron_heating = Grid::Zero(m_xdim,m_ydim);
    m_cumulative_ion_heating = m_cumulative_electron_heating;
    m_cumulative_joule_heating = m_cumulative_electron_heating;
  }
  
  // overwrite the init.state file if new simulation with overwrite flag
  if(!continue_mode && m_overwrite_init){
    std::cout << "Writing out init.state...\n";
    writeStateFile("init");
  }
  // initialize container for data to write to mhd.out
  initOutputContainer();
  // initialize the mhd.out file
  if(!continue_mode && m_write_interval > 0){
    outputPreamble();
    storeGrids();
    writeToOutFile();
  }
  
}

// Initialoize the size of m_data_to_write and the capacity of each std::string element
void PlasmaDomain::initOutputContainer()
{
  // initialize grid with copies of pi, to simulate largest precision necessary for allocating space
  Grid pi_grid(m_xdim,m_ydim,PI*1e100); // this represents the largest number we expect to write to mhd.out
  std::string row = pi_grid.format(',','\n');
  // count number of lines to record per iteration
    // one line for time
    // one line for each grid name written from equation set or module
    // m_xdim lines for grid written from equation set or module
  int num_grids_to_record{0};
  for(int i=0; i<m_eqs->num_variables(); i++){
    if(m_eqs->getOutputFlag(i)) num_grids_to_record++;
  }
  std::vector<std::string> module_varnames;
  std::vector<Grid> module_data;
  m_module_handler.getFileOutputData(module_varnames,module_data);
  for(int i=0; i<module_varnames.size(); i++){
    num_grids_to_record++;
  }
  int num_lines_per_iter = 1+2*num_grids_to_record;
  if(m_multispecies_mode) num_lines_per_iter += 6; //two each for cumulative_electron_heating and cumulative_ion_heating and cumulative_joule_heating
  // initialize size of m_data_to_write and the capacity of each element
  if (m_write_interval > 0){
    m_data_to_write.resize(m_write_interval*num_lines_per_iter);
    for (auto& elem : m_data_to_write) elem.reserve(row.size());
  }
}

int PlasmaDomain::gridname2index(const std::string& name) const
{
  auto it = std::find(m_gridnames.begin(),m_gridnames.end(),name);
  if (it==m_gridnames.end()){
    std::cerr << "<" << name << "> is not a valid grid name." << std::endl;
    assert(false);
  }
  return it - m_gridnames.begin();
}

bool PlasmaDomain::is_grid(const std::string& name) const
{
  auto it = std::find(m_gridnames.begin(),m_gridnames.end(),name);
  return it != m_gridnames.end();
}

Grid& PlasmaDomain::grid(int index)
{
  assert(index>=0 && index<m_gridnames.size());
  return m_grids[index];
}

Grid& PlasmaDomain::grid(const std::string& name)
{
  int index = gridname2index(name);
  return m_grids[index];
}

//Compute lower and upper x- and y- indicies for differential operations
//from boundary conditions, to exclude ghost zones. Results stored in m_xl, m_xu, m_yl, m_yu
//Also generates m_ghost_zone_mask Grid.
//Also determines m_xl_dt, m_xu_dt, m_yl_dt, m_yu_dt, i.e. bounds for timestep computation
void PlasmaDomain::computeIterationBounds()
{
  assert(m_xdim > 2*N_GHOST && m_ydim > 2*N_GHOST && "Grid too small for ghost zones");
  if(x_bound_1 == BoundaryCondition::Open || x_bound_1 == BoundaryCondition::Reflect || x_bound_1 == BoundaryCondition::Fixed || x_bound_1 == BoundaryCondition::OpenMoC || x_bound_1 == BoundaryCondition::OpenUCNP) m_xl = N_GHOST;
  else{ assert(x_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xl = 0; }
  if(x_bound_2 == BoundaryCondition::Open || x_bound_2 == BoundaryCondition::Reflect || x_bound_2 == BoundaryCondition::Fixed || x_bound_2 == BoundaryCondition::OpenMoC || x_bound_2 == BoundaryCondition::OpenUCNP) m_xu = m_xdim - N_GHOST - 1;
  else{ assert(x_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_xu = m_xdim - 1; }
  if(y_bound_1 == BoundaryCondition::Open || y_bound_1 == BoundaryCondition::Reflect || y_bound_1 == BoundaryCondition::Fixed || y_bound_1 == BoundaryCondition::OpenMoC || y_bound_1 == BoundaryCondition::OpenUCNP) m_yl = N_GHOST;
  else{ assert(y_bound_1 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yl = 0; }
  if(y_bound_2 == BoundaryCondition::Open || y_bound_2 == BoundaryCondition::Reflect || y_bound_2 == BoundaryCondition::Fixed || y_bound_2 == BoundaryCondition::OpenMoC || y_bound_2 == BoundaryCondition::OpenUCNP) m_yu = m_ydim - N_GHOST - 1;
  else{ assert(y_bound_2 == BoundaryCondition::Periodic && "Boundary cond'n must be defined"); m_yu = m_ydim - 1; }
  m_ghost_zone_mask = Grid(m_xdim,m_ydim,0.0);
  for(int i = m_xl; i <= m_xu; i++){
    for(int j = m_yl; j <= m_yu; j++){
      m_ghost_zone_mask(i,j) = 1.0;
    }
  }
  m_xl_dt = m_xl, m_yl_dt = m_yl, m_xu_dt = m_xu, m_yu_dt = m_yu;
  if(x_bound_1==BoundaryCondition::OpenMoC) m_xl_dt -= N_GHOST;
  if(x_bound_2==BoundaryCondition::OpenMoC) m_xu_dt += N_GHOST;
  if(y_bound_1==BoundaryCondition::OpenMoC) m_yl_dt -= N_GHOST;
  if(y_bound_2==BoundaryCondition::OpenMoC) m_yu_dt += N_GHOST;

}

//Takes not-necessarily-uniform cell sizes d_x and d_y and converts to cell center positions,
  //returned as the first (x position) and second (y position) indicies of a vector of Grids
  //x_origin and y_origin specify location of the origin (0,0) in the domain:
  //"lower" is the default, "center" and "upper" are also options for each
  //For example, x_origin = "center" and y_origin = "center" places the origin
  //at the center of the domain.
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
      double diff_pos = pos(i_next,j_next) - pos(i,j);
      double diff_d = 0.5*d(i,j) + 0.5*d(i_next,j_next);
      valid = valid && ( std::abs( 1.0 - diff_pos / diff_d ) < tolerance);
    }
  }
  return valid;
}

bool PlasmaDomain::allInternalGridsInitialized()
{
  for(const Grid& g : m_grids){
    if(g.size() == 1) return false;
  }
  return true;
}
