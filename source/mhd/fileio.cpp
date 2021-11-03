//fileio.cpp
//PlasmaDomain functionality relating to file I/O

#include "plasmadomain.hpp"

//Must match ordering of BoundaryCondition enum in plasmadomain.hpp
const std::vector<std::string> PlasmaDomain::m_boundary_condition_names = {
  "periodic", "open", "fixed", "reflect"
};

//Corresponding variable names for file I/O
//Must match the ordering of the Variable enum defined in plasmadomain.hpp
const std::vector<std::string> PlasmaDomain::m_var_names = {
  "d_x","d_y","pos_x","pos_y","rho","temp","mom_x","mom_y","mom_z","b_x","b_y","b_z","grav_x","grav_y",
  "press","thermal_energy","kinetic_energy","rad","dt","dt_thermal","dt_rad","v_x","v_y","v_z","n",
  "b_magnitude","b_hat_x","b_hat_y",
  "mag_press","lorentz_force_x","lorentz_force_y","mag_pxx","mag_pyy","mag_pzz","mag_pxy","mag_pxz","mag_pyz"
};

//For purposes of reading in .config files
//Must match ordering of Config enum below
const std::vector<std::string> PlasmaDomain::m_config_names = {
  "x_bound_1","x_bound_2","y_bound_1","y_bound_2","radiative_losses","ambient_heating",
  "thermal_conduction","flux_saturation","temp_chromosphere","radiation_ramp","heating_rate",
  "epsilon","epsilon_thermal","epsilon_rad","epsilon_viscous","dt_thermal_min","rho_min",
  "temp_min","thermal_energy_min","max_iterations","iter_output_interval","time_output_interval",
  "output_flags","xdim","ydim","open_boundary_strength","std_out_interval","safe_state_mode",
  "open_boundary_decay_base", "x_origin", "y_origin"
};

enum class Config {
  x_bound_1, x_bound_2, y_bound_1, y_bound_2, radiative_losses, ambient_heating,
  thermal_conduction, flux_saturation, temp_chromosphere, radiation_ramp, heating_rate,
  epsilon, epsilon_thermal, epsilon_rad, epsilon_viscous, dt_thermal_min, rho_min,
  temp_min, thermal_energy_min, max_iterations, iter_output_interval, time_output_interval,
  output_flags, xdim, ydim, open_boundary_strength, std_out_interval, safe_state_mode,
  open_boundary_decay_base, x_origin, y_origin
};

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
  assert(line == "xdim,ydim");
  std::getline(in_file, line);
  std::istringstream ss_dim(line);
  std::getline(ss_dim,element,',');
  m_xdim = std::stoi(element);
  std::getline(ss_dim,element,',');
  m_ydim = std::stoi(element);
  for(int var=0; var<num_variables; var++) m_grids[var] = Grid(m_xdim,m_ydim,0.0);

  std::getline(in_file, line);
  assert(line == "ion_mass");
  std::getline(in_file, line);
  m_ion_mass = std::stod(line);

  std::getline(in_file, line);
  assert(line == "adiabatic_index");
  std::getline(in_file, line);
  m_adiabatic_index = std::stod(line);

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
  computeIterationBounds();
  computeConstantTerms();
  recomputeDerivedVariables();
}

//Read in simulation configuration from .config file
//Allows to change configuration without recompiling
void PlasmaDomain::readConfigFile(const char* config_filename)
{
  std::ifstream in_file(config_filename);
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
  m_out_file << "xdim,ydim" << std::endl;
  m_out_file << m_xdim << "," << m_ydim << std::endl;
  m_out_file << "pos_x" << std::endl;
  m_out_file << m_grids[pos_x].format(',',';');
  m_out_file << "pos_y" << std::endl;
  m_out_file << m_grids[pos_y].format(',',';');
  m_out_file << "b_x" << std::endl;
  m_out_file << m_grids[b_x].format(',',';');
  m_out_file << "b_y" << std::endl;
  m_out_file << m_grids[b_y].format(',',';');
  m_out_file << "b_z" << std::endl;
  m_out_file << m_grids[b_z].format(',',';');
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
//All variables included in PlasmaDomain::StateVars are written out here
void PlasmaDomain::writeStateFile(std::string filename,int precision) const
{
  std::ofstream state_file;
  if (filename.compare("mhd") == 0){
    if(safe_state_mode) state_file.open(m_out_directory/(filename+std::to_string(m_iter%2)+".state"));
    else state_file.open(m_out_directory/(filename+".state"));
  }
  else state_file.open(m_out_directory/(filename+".state"));
  state_file << "xdim,ydim\n";
  state_file << m_xdim << "," << m_ydim << std::endl;
  state_file << "ion_mass\n";
  state_file << m_ion_mass << std::endl;
  state_file << "adiabatic_index\n";
  state_file << m_adiabatic_index << std::endl;
  state_file << "t=" << m_time << std::endl;
  for(int i=state_var_start; i<state_var_end; i++){
    state_file << m_var_names[i] << std::endl;
    state_file << m_grids[i].format(',',';',precision);
  }
  state_file.close();
}

//Gets rid of the older of the two state files at the end of the sim run
void PlasmaDomain::cleanUpStateFiles() const
{
  std::filesystem::rename(m_out_directory/("mhd"+std::to_string((m_iter)%2)+".state"),
          m_out_directory/("mhd.state"));
  std::filesystem::remove(m_out_directory/("mhd"+std::to_string((m_iter-1)%2)+".state"));
}

PlasmaDomain::BoundaryCondition PlasmaDomain::stringToBoundaryCondition(const std::string str) const
{
  auto it = std::find(m_boundary_condition_names.begin(),m_boundary_condition_names.end(),str);
  assert(it != m_boundary_condition_names.end()); //Ensure that given name is valid boundary condition
  auto index = std::distance(m_boundary_condition_names.begin(),it);
  return static_cast<BoundaryCondition>(index);
}

void PlasmaDomain::printUpdate(double min_dt, int subcycles_thermal, int subcycles_rad) const
{
  // printf("Iter: %i",m_iter);
  // if(max_iterations > 0) printf("/%i",max_iterations);
  // printf("|t: %f",m_time);
  // if(max_time > 0.0) printf("/%f",max_time);
  // printf("|dt: %f",min_dt);
  // if(thermal_conduction) printf("|Cond. Subcyc: %i",subcycles_thermal);
  // if(radiative_losses) printf("|Rad. Subcyc: %i",subcycles_rad);
  // printf("\n");
  std::cout << "Iter: " << m_iter;
  if(max_iterations > 0) std::cout << "/" << max_iterations;
  std::cout << "|t: " << m_time;
  if(max_time > 0.0) std::cout << "/" << max_time;
  std::cout << "|dt: " << min_dt;
  if(thermal_conduction) std::cout << "|Cond. Subcyc: " << subcycles_thermal;
  if(radiative_losses) std::cout << "|Rad. Subcyc: " << subcycles_rad;
  std::cout << std::endl;
}

void PlasmaDomain::handleSingleConfig(int setting_index, std::string rhs)
{
  switch (static_cast<int>(setting_index)) {
  case static_cast<int>(Config::x_bound_1): x_bound_1 = stringToBoundaryCondition(rhs); break;
  case static_cast<int>(Config::x_bound_2): x_bound_2 = stringToBoundaryCondition(rhs); break;
  case static_cast<int>(Config::y_bound_1): y_bound_1 = stringToBoundaryCondition(rhs); break;
  case static_cast<int>(Config::y_bound_2): y_bound_2 = stringToBoundaryCondition(rhs); break;
  case static_cast<int>(Config::radiative_losses): radiative_losses = (rhs == "true"); break;
  case static_cast<int>(Config::ambient_heating): ambient_heating = (rhs == "true"); break;
  case static_cast<int>(Config::thermal_conduction): thermal_conduction = (rhs == "true"); break;
  case static_cast<int>(Config::flux_saturation): flux_saturation = (rhs == "true"); break;
  case static_cast<int>(Config::temp_chromosphere): temp_chromosphere = std::stod(rhs); break;
  case static_cast<int>(Config::radiation_ramp): radiation_ramp = std::stod(rhs); break;
  case static_cast<int>(Config::heating_rate): heating_rate = std::stod(rhs); break;
  case static_cast<int>(Config::epsilon): epsilon = std::stod(rhs); break;
  case static_cast<int>(Config::epsilon_thermal): epsilon_thermal = std::stod(rhs); break;
  case static_cast<int>(Config::epsilon_rad): epsilon_rad = std::stod(rhs); break;
  case static_cast<int>(Config::epsilon_viscous): epsilon_viscous = std::stod(rhs); break;
  case static_cast<int>(Config::dt_thermal_min): dt_thermal_min = std::stod(rhs); break;
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
  case static_cast<int>(Config::std_out_interval): std_out_interval = std::stoi(rhs); break;
  case static_cast<int>(Config::open_boundary_decay_base): open_boundary_decay_base = std::stod(rhs); break;
  case static_cast<int>(Config::x_origin): x_origin = rhs; break;
  case static_cast<int>(Config::y_origin): y_origin = rhs; break;
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