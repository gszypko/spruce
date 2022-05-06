//fileio.cpp
//PlasmaDomain functionality relating to file I/O

#include "plasmadomain.hpp"

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
    if(line[0] == '#') comment_lines.push_back(line);
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
  for(int var=0; var<num_variables; var++) m_grids.push_back(Grid(m_xdim,m_ydim,0.0));

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

  bool first_loop = true;
  while(getCleanedLine(in_file, line)){
    if(first_loop){
      if(line.find("duration")!=std::string::npos){
        std::istringstream ss_line(line);
        std::string el;
        std::getline(ss_line,el,'=');
        assert(el == "duration" && "Duration specification must be formatted as duration=[duration]");
        std::getline(ss_line,el);
        clearWhitespace(el);
        assert( (std::isdigit(el[0]) || el[0] == '.' || el[0] == '-') && "Encountered non-numerical value for duration from state file");
        m_duration = std::stod(el);
        first_loop = false;
        continue;
      }
      first_loop = false;
    }
    //Ensure that variable in file is valid
    auto it = std::find(m_var_names.begin(),m_var_names.end(),line);
    if(it == m_var_names.end()){
      std::cerr << "Variable name " << line << " in state file not recognized\n";
      abort();
    }
    auto index = std::distance(m_var_names.begin(),it);
    Grid curr_grid(m_xdim,m_ydim);
    //Read in Grid corresponding to variable
    int j;
    std::string row; std::string el;
    for(int i=0; i<m_xdim; i++){
      j=0;
      getCleanedLine(in_file,row);
      std::istringstream ss_row(row);
      assert((std::isdigit(ss_row.peek()) || (ss_row.peek()=='-') || (ss_row.peek()=='.')) && "Encountered non-numerical row in .state file sooner than expected");
      while(std::getline(ss_row,el,',')){
        assert(j < m_ydim && "Row in .state file is too long (greater than ydim)");
        curr_grid(i,j) = std::stod(el);
        j++;
      }
    }
    assert(!std::isdigit(in_file.peek()) && !(in_file.peek()=='-') && "Encountered more rows in a .state file grid than expected");
    m_grids[index] = curr_grid;
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
  //Read through config file
  while(std::getline(in_file, line)){
    clearWhitespace(line);
    if(line.empty() || line[0] == '#') continue; //skip comment and empty lines
    std::istringstream ss_line(line);
    std::string lhs, rhs;
    std::getline(ss_line,lhs,'=');
    std::getline(ss_line,rhs,'#');
    if(m_module_handler.isModuleName(lhs)){
      m_module_handler.instantiateModule(lhs,in_file,(rhs == "true"));
      continue;
    }
    std::vector<std::string> rhs_vec = splitString(rhs,',');
    auto it = std::find(m_config_names.begin(),m_config_names.end(),lhs);
    assert(it != m_config_names.end());
    auto index = std::distance(m_config_names.begin(),it);
    if(rhs_vec.size() == 1) handleSingleConfig(index,rhs);
    else if(rhs_vec.size() > 1) handleConfigList(index,rhs_vec,rhs_lists,list_vars,num_combinations);
  }
}

void PlasmaDomain::outputPreamble()
{
  for(std::string comment : comment_lines){
    assert(comment[0] == '#');
    m_out_file << comment << std::endl;
  }
  m_out_file << "xdim,ydim" << std::endl;
  m_out_file << m_xdim << "," << m_ydim << std::endl;
  m_out_file << "pos_x" << std::endl;
  m_out_file << m_grids[pos_x].format(',',';');
  m_out_file << "pos_y" << std::endl;
  m_out_file << m_grids[pos_y].format(',',';');
  m_out_file << "be_x" << std::endl;
  m_out_file << m_grids[be_x].format(',',';');
  m_out_file << "be_y" << std::endl;
  m_out_file << m_grids[be_y].format(',',';');
}

void PlasmaDomain::outputCurrentState()
{
  m_out_file << "t=" << m_time << std::endl;
  for(int i=0; i<num_variables; i++){
    if(m_output_flags[i]) writeGridToOutput(m_grids[i],m_var_names[i]);
  }
  std::vector<std::string> module_varnames;
  std::vector<Grid> module_data;
  m_module_handler.getFileOutputData(module_varnames,module_data);
  for(int i=0; i<module_varnames.size(); i++){
    writeGridToOutput(module_data[i],module_varnames[i]);
  }
}

void PlasmaDomain::writeGridToOutput(const Grid& grid, std::string var_name)
{
  m_out_file << var_name << std::endl;
  m_out_file << grid.format(',',';');
}

//Output current state into state file (toggles between overwriting two different files
//so that the most recent state is still written in event of a crash)
//All variables included in PlasmaDomain::StateVars are written out here
//Resulting file is named [filename_stem].state
void PlasmaDomain::writeStateFile(std::string filename_stem,int precision) const
{
  std::ofstream state_file;
  if (filename_stem.compare("mhd") == 0){
    if(safe_state_mode) state_file.open(m_out_directory/fs::path(filename_stem+std::to_string((m_iter/safe_state_interval)%2)+".state"));
    else state_file.open(m_out_directory/fs::path(filename_stem+".state"));
  }
  else state_file.open(m_out_directory/fs::path(filename_stem+".state"));
  for(std::string comment : comment_lines){
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
  for(int i : state_vars){
    state_file << m_var_names[i] << std::endl;
    state_file << m_grids[i].format(',','\n',precision);
  }
  state_file.close();
}

//Gets rid of the older of the two state files, named mhd0.state and mhd1.state, at the end of the sim run
//The remaining one is renamed to [filename_stem].state
void PlasmaDomain::cleanUpStateFiles(std::string filename_stem) const
{
  fs::rename(m_out_directory/("mhd"+std::to_string((m_iter/safe_state_interval)%2)+".state"),
          m_out_directory/(filename_stem+".state"));
  fs::remove(m_out_directory/("mhd"+std::to_string((m_iter/safe_state_interval - 1)%2)+".state"));
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
  if(max_time > 0.0) std::cout << "/" << max_time;
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
  case static_cast<int>(Config::rho_min): rho_min = std::stod(rhs); break;
  case static_cast<int>(Config::temp_min): temp_min = std::stod(rhs); break;
  case static_cast<int>(Config::thermal_energy_min): thermal_energy_min = std::stod(rhs); break;
  case static_cast<int>(Config::max_iterations): max_iterations = std::stoi(rhs); break;
  case static_cast<int>(Config::iter_output_interval): iter_output_interval = std::stoi(rhs); break;
  case static_cast<int>(Config::time_output_interval): time_output_interval = std::stod(rhs); break;
  case static_cast<int>(Config::output_flags): setOutputFlag(rhs,true); break;
  case static_cast<int>(Config::xdim): m_xdim = std::stoi(rhs); break;
  case static_cast<int>(Config::ydim): m_ydim = std::stoi(rhs); break;
  case static_cast<int>(Config::open_boundary_strength): open_boundary_strength = std::stod(rhs); break;
  case static_cast<int>(Config::safe_state_mode): safe_state_mode = (rhs == "true"); break;
  case static_cast<int>(Config::safe_state_interval): safe_state_interval = std::stoi(rhs); break;
  case static_cast<int>(Config::std_out_interval): std_out_interval = std::stoi(rhs); break;
  case static_cast<int>(Config::open_boundary_decay_base): open_boundary_decay_base = std::stod(rhs); break;
  // case static_cast<int>(Config::x_origin): x_origin = rhs; break;
  // case static_cast<int>(Config::y_origin): y_origin = rhs; break;
  case static_cast<int>(Config::time_integrator): time_integrator = stringToTimeIntegrator(rhs); break;
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

void PlasmaDomain::setOutputFlag(std::string var_name, bool new_flag)
{
  auto it = std::find(m_var_names.begin(),m_var_names.end(),var_name);
  assert(it != m_var_names.end());
  auto index = std::distance(m_var_names.begin(),it);
  m_output_flags[index] = new_flag;
}

void PlasmaDomain::setOutputFlags(const std::vector<std::string> var_names, bool new_flag)
{
  for(std::string var_name : var_names) setOutputFlag(var_name,new_flag);
}