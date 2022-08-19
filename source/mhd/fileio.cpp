//fileio.cpp
//PlasmaDomain functionality relating to file I/O

#include "plasmadomain.hpp"
#include "utils.hpp"
#include <iostream>

//Read in variables from .state file
//This function will abort execution if an invalid variable name is encountered
//Does not check that the variables read from the file are a complete
//description of the plasma, nor that none of them are contradictory
//Any lines at the top beginning with '#' will be stored as comment lines
//and written out at the top of all .state and .out files for documentation purposes
void PlasmaDomain::readStateFile(const fs::path &state_file, bool continue_mode)
{  
  std::ifstream in_file(state_file.string());

  //Read in grid dimensions
  std::string line, element;
  std::getline(in_file,line);
  while(line.empty() || line[0] == '#'){
    if(line[0] == '#') m_comment_lines.push_back(line);
    std::getline(in_file,line);
  }
  clearWhitespace(line);
  assert(line == "xdim,ydim");
  getCleanedLine(in_file, line);
  std::istringstream ss_dim(line);
  std::getline(ss_dim,element,',');
  m_xdim = std::stoi(element);
  std::getline(ss_dim,element,',');
  m_ydim = std::stoi(element);

  getCleanedLine(in_file, line);
  assert(line == "ion_mass");
  getCleanedLine(in_file, line);
  m_ion_mass = std::stod(line);

  getCleanedLine(in_file, line);
  assert(line == "adiabatic_index");
  getCleanedLine(in_file, line);
  m_adiabatic_index = std::stod(line);

  //Read in time of state file
  getCleanedLine(in_file,line);
  std::istringstream ss_time(line);
  std::getline(ss_time,element,'=');
  assert(element == "t");
  std::getline(ss_time,element);
  if(continue_mode) m_time = std::stod(element);
  else m_time = 0.0;

  while(getCleanedLine(in_file, line)){
    std::string var_name = line;
    Grid curr_grid(m_xdim,m_ydim);
    //Read in Grid corresponding to variable
    int j;
    std::string row; std::string el;
    for(int i=0; i<m_xdim; i++){
      j=0;
      getCleanedLine(in_file,row);
      std::istringstream ss_row(row);
      assert((std::isdigit(ss_row.peek()) || (ss_row.peek()=='-') || (ss_row.peek()=='.') || (ss_row.peek()=='n') || (ss_row.peek()=='i'))
              && "Encountered non-numerical row in .state file sooner than expected");
      while(std::getline(ss_row,el,',')){
        assert(j < m_ydim && "Row in .state file is too long (greater than ydim)");
        curr_grid(i,j) = std::stod(el);
        j++;
      }
    }
    assert(!std::isdigit(in_file.peek()) && !(in_file.peek()=='-') && "Encountered more rows in a .state file grid than expected");
    auto it = std::find(m_gridnames.begin(),m_gridnames.end(),var_name);
    if(it != m_gridnames.end()){
      auto index = std::distance(m_gridnames.begin(),it);
      m_grids[index] = curr_grid;
    }
    else m_eqs->grid(var_name) = curr_grid;
  }
  in_file.close();
}

//Read in simulation configuration from .config file
//Allows to change configuration without recompiling
void PlasmaDomain::readConfigFile(const fs::path &config_file)
{
  std::ifstream in_file(config_file.string());
  std::string line;
  std::vector<std::vector<std::string> > rhs_lists;
  std::vector<int> list_vars;
  int num_combinations = 1;
  // read through config file
  while (std::getline(in_file, line)){
    // obtain the lhs and rhs quantities for the current line
    clearWhitespace(line);
    if (line.empty() || line[0] == '#') continue; // skip comment and empty lines
    std::istringstream ss_line(line);
    std::string lhs, rhs;
    std::getline(ss_line,lhs,'=');
    std::getline(ss_line,rhs,'#');
    // check if lhs corresponds to the equation set and, if so, instantiate that equation set
    if (EquationSet::isEquationSetName(lhs)){
      m_eqs = EquationSet::instantiateWithConfig((*this),in_file,lhs,(rhs == "true"));
      continue;
    }
    // check if lhs corresponds to a module and, if so, instantiate the module
    if (m_module_handler.isModuleName(lhs)){
      assert(m_eqs.get() != nullptr && "<m_eqs> must be instantiated before instantiating modules");
      m_module_handler.instantiateModule(lhs,in_file,(rhs == "true"));
      continue;
    }
    // if config does not correspond to equation set or module, handle as PlasmaDomain config
    std::vector<std::string> rhs_vec = splitString(rhs,',');
    auto it = std::find(m_config_names.begin(),m_config_names.end(),lhs);
    if (it == m_config_names.end()){
      std::cerr << lhs << " is not a valid config name." << std::endl;
      assert(false);
    }
    auto index = std::distance(m_config_names.begin(),it);
    if(rhs_vec.size() == 1) handleSingleConfig(index,rhs);
    else if(rhs_vec.size() > 1) handleConfigList(index,rhs_vec,rhs_lists,list_vars,num_combinations);
  }
}

// initializes mhd.out with the comment lines and internal grid information
void PlasmaDomain::outputPreamble()
{
  std::ofstream out_file(m_out_directory/m_out_filename);
  for (std::string comment : m_comment_lines){
    assert(comment[0] == '#');
    out_file << comment << std::endl;
  }
  out_file << "xdim,ydim" << std::endl;
  out_file << m_xdim << "," << m_ydim << std::endl;
  for (int v : {pos_x,pos_y,be_x,be_y}){
    out_file << m_gridnames[v] << std::endl;
    out_file << m_grids[v].format(',','\n');
  }
  out_file.close();
}

// Stores the current time and grids from the equation set and modules within m_data_to_write
  // Each element of m_data_to_write that is updated by this function should be empty before it
  // is updated. 
void PlasmaDomain::storeGrids()
{
  // write time
  assert(m_lines_recorded<m_data_to_write.size() && m_data_to_write[m_lines_recorded].empty());
  m_data_to_write[m_lines_recorded] = "t=" + num2str(m_time) + '\n';
  m_lines_recorded++;

  // write grids from equation set that have output flag
  for (int i=0; i<m_eqs->num_variables(); i++){
    if (m_eqs->getOutputFlag(i)){
      assert(m_lines_recorded<m_data_to_write.size() && m_data_to_write[m_lines_recorded].empty());
      m_data_to_write[m_lines_recorded] = m_eqs->nameFromIndex(i) + '\n';
      m_lines_recorded++;
      assert(m_lines_recorded<m_data_to_write.size() && m_data_to_write[m_lines_recorded].empty());
      m_data_to_write[m_lines_recorded] = m_eqs->grid(i).format(',','\n');
      m_lines_recorded++;
    }
  }

  // write grids from modules
  std::vector<std::string> module_varnames;
  std::vector<Grid> module_data;
  m_module_handler.getFileOutputData(module_varnames,module_data);
  for(int i=0; i<module_varnames.size(); i++){
    assert(m_lines_recorded<m_data_to_write.size() && m_data_to_write[m_lines_recorded].empty());
    m_data_to_write[m_lines_recorded] = module_varnames[i] + '\n';
    m_lines_recorded++;
    assert(m_lines_recorded<m_data_to_write.size() && m_data_to_write[m_lines_recorded].empty());
    m_data_to_write[m_lines_recorded] = module_data[i].format(',','\n');
    m_lines_recorded++;
  }
  // increment store counter
  m_store_counter++;
}

// Writes the information stored within m_data_to_write to mhd.out
  // Each element of m_data_to_write is cleared after it is written. This preserves the capacity
  // of the std::string element while ensuring that the same data is not written twice.
void PlasmaDomain::writeToOutFile()
{
  std::ofstream out_file(m_out_directory/m_out_filename,std::ofstream::app);
  for (int line = 0; line<m_lines_recorded; line++){
    assert(!m_data_to_write[line].empty());
    out_file << m_data_to_write[line];
    m_data_to_write[line].clear();
  }
  m_lines_recorded = 0;
  m_store_counter = 0;
}

//Output current state into state file (toggles between overwriting two different files
//so that the most recent state is still written in event of a crash)
//All variables included in PlasmaDomain::StateVars are written out here
//Resulting file is named [filename_stem].state
void PlasmaDomain::writeStateFile(std::string filename_stem,int precision) const
{
  // determine name of state file
  fs::path filename;
  if (filename_stem == "mhd"){
    filename = "mhd" + num2str(m_state_identifier) + ".state";
  }
  else filename = filename_stem + ".state";

  // open stream to state file
  std::ofstream state_file;
  state_file.open(m_out_directory/filename);
  
  // write to state file
  for(std::string comment : m_comment_lines){
    assert(comment[0] == '#');
    state_file << comment << std::endl;
  }
  state_file << "xdim,ydim\n";
  state_file << m_xdim << "," << m_ydim << std::endl;
  state_file << "ion_mass\n";
  state_file << m_ion_mass << std::endl;
  state_file << "adiabatic_index\n";
  state_file << m_adiabatic_index << std::endl;
  state_file << "t=" << m_time << std::endl;
  for(int i=0; i<m_gridnames.size(); i++){
    state_file << m_gridnames[i] << std::endl;
    state_file << m_grids[i].format(',','\n',precision);
  }
  for(int i : m_eqs->state_variables()){
    state_file << m_eqs->nameFromIndex(i) << std::endl;
    state_file << m_eqs->grid(i).format(',','\n',precision);
  }
  state_file.close();
}

void PlasmaDomain::updateStateIdentifier()
{
  if (m_state_identifier == 1) m_state_identifier = 2;
  else if (m_state_identifier == 2) m_state_identifier = 1;
  else{
    std::cerr << "<m_state_identifier> out of range: must be 1 or 2" << std::endl;
    assert(false);
  }
}

PlasmaDomain::BoundaryCondition PlasmaDomain::stringToBoundaryCondition(const std::string str) const
{
  auto it = std::find(m_boundary_condition_names.begin(),m_boundary_condition_names.end(),str);
  assert(it != m_boundary_condition_names.end()); //Ensure that given name is valid boundary condition
  auto index = std::distance(m_boundary_condition_names.begin(),it);
  return static_cast<BoundaryCondition>(index);
}

PlasmaDomain::TimeIntegrator PlasmaDomain::stringToTimeIntegrator(const std::string str) const
{
  auto it = std::find(m_time_integrator_names.begin(),m_time_integrator_names.end(),str);
  assert(it != m_time_integrator_names.end()); //Ensure that given name is valid
  auto index = std::distance(m_time_integrator_names.begin(),it);
  return static_cast<TimeIntegrator>(index);
}

void PlasmaDomain::printUpdate(double dt) const
{
  std::cout << "Iter: " << m_iter;
  if(max_iterations > 0) std::cout << "/" << max_iterations;
  std::cout << "|t: " << m_time;
  if(m_max_time > 0.0) std::cout << "/" << m_max_time;
  std::cout << "|dt: " << dt;
  std::vector<std::string> module_messages = m_module_handler.getCommandLineMessages();
  for(int i=0; i<module_messages.size(); i++){
    std::cout << "|" << module_messages[i];
  }
  std::cout << std::endl;
}

void PlasmaDomain::handleSingleConfig(int setting_index, std::string rhs)
{
  switch (static_cast<int>(setting_index)) {
    case static_cast<int>(Config::x_bound_1): x_bound_1 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Config::x_bound_2): x_bound_2 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Config::y_bound_1): y_bound_1 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Config::y_bound_2): y_bound_2 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Config::epsilon): epsilon = std::stod(rhs); break;
    case static_cast<int>(Config::epsilon_viscous): epsilon_viscous = std::stod(rhs); break;
    case static_cast<int>(Config::density_min): density_min = std::stod(rhs); break;
    case static_cast<int>(Config::temp_min): temp_min = std::stod(rhs); break;
    case static_cast<int>(Config::thermal_energy_min): thermal_energy_min = std::stod(rhs); break;
    case static_cast<int>(Config::max_iterations): max_iterations = std::stoi(rhs); break;
    case static_cast<int>(Config::iter_output_interval): m_iter_output_interval = std::stoi(rhs); break;
    case static_cast<int>(Config::time_output_interval): m_time_output_interval = std::stod(rhs); break;
    case static_cast<int>(Config::output_flags): assert(m_eqs.get() != nullptr && "equation_set must be defined in config file before output_flags"); m_eqs->setOutputFlag(rhs,true); break;
    case static_cast<int>(Config::xdim): m_xdim = std::stoi(rhs); break;
    case static_cast<int>(Config::ydim): m_ydim = std::stoi(rhs); break;
    case static_cast<int>(Config::open_boundary_strength): open_boundary_strength = std::stod(rhs); break;
    case static_cast<int>(Config::write_interval): m_write_interval = std::stoi(rhs); break;
    case static_cast<int>(Config::std_out_interval): m_std_out_interval = std::stoi(rhs); break;
    case static_cast<int>(Config::open_boundary_decay_base): open_boundary_decay_base = std::stod(rhs); break;
    case static_cast<int>(Config::time_integrator): m_time_integrator = stringToTimeIntegrator(rhs); break;
    case static_cast<int>(Config::duration): m_duration = std::stod(rhs); break;
    case static_cast<int>(Config::sg_opt): m_sg_opt = rhs; break;
    default: break;
  }
}

void PlasmaDomain::handleConfigList(int setting_index, std::vector<std::string> rhs_vec,
                                     std::vector<std::vector<std::string> > &rhs_lists,
                                     std::vector<int> &list_vars, int &num_combinations)
{
  if(setting_index == static_cast<int>(Config::output_flags)) setOutputFlags(rhs_vec,true);
  // else if(setting_index == static_cast<int>(Config::state_flags)) setStateFlags(rhs_vec,true);
  else assert(false && "Only output_flags support list specifier");
}

void PlasmaDomain::setOutputFlags(const std::vector<std::string> var_names, bool new_flag)
{
  assert(m_eqs.get() != nullptr && "equation_set must be defined in config file before output_flags");
  for(std::string var_name : var_names) m_eqs->setOutputFlag(var_name,new_flag);
}