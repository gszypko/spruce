//fileio.cpp
//PlasmaDomain functionality relating to file I/O

#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "plasmadomain.hpp"

//Must match ordering of BoundaryCondition enum in plasmadomain.hpp
const std::vector<std::string> PlasmaDomain::m_boundary_condition_names = {
  "periodic","wall","open"
};

//Corresponding variable names for file I/O
//Must match the ordering of the Variable enum defined in plasmadomain.hpp
const std::vector<std::string> PlasmaDomain::m_var_names = {
  "rho","temp","mom_x","mom_y","b_x","b_y","b_z",
  "press","energy","rad","dt","dt_thermal","dt_rad","v_x","v_y",
  "grav_x","grav_y","b_magnitude","b_hat_x","b_hat_y",
  "mag_press","mag_pxx","mag_pyy","mag_pzz","mag_pxy","mag_pxz","mag_pyz"
};

//For purposes of reading in .settings files
//Must match ordering of Settings enum below
const std::vector<std::string> PlasmaDomain::m_setting_names = {
  "x_bound_1","x_bound_2","y_bound_1","y_bound_2","radiative_losses","ambient_heating",
  "thermal_conduction","flux_saturation","temp_chromosphere","radiation_ramp","heating_rate",
  "b_0","epsilon","epsilon_thermal","epsilon_rad","epsilon_viscous","dt_thermal_min","rho_min",
  "temp_min","thermal_energy_min","max_iterations","max_time","iter_output_interval","time_output_interval",
  "output_flags","state_flags","xdim","ydim","dx","dy"
};

enum class Setting {
  x_bound_1, x_bound_2, y_bound_1, y_bound_2, radiative_losses, ambient_heating,
  thermal_conduction, flux_saturation, temp_chromosphere, radiation_ramp, heating_rate,
  b_0, epsilon, epsilon_thermal, epsilon_rad, epsilon_viscous, dt_thermal_min, rho_min,
  temp_min, thermal_energy_min, max_iterations, max_time, iter_output_interval, time_output_interval,
  output_flags, state_flags, xdim, ydim, dx, dy
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
  std::istringstream ss_dim(line);
  std::getline(ss_dim,element,',');
  m_xdim = std::stoi(element);
  std::getline(ss_dim,element,',');
  m_ydim = std::stoi(element);

  //Read in physical dimensions
  std::getline(in_file,line);
  std::istringstream ss_dxdy(line);
  std::getline(ss_dxdy,element,',');
  m_dx = std::stod(element);
  std::getline(ss_dxdy,element,',');
  m_dy = std::stod(element);

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
  computeMagneticTerms();
  recomputeDerivedVariables();
}

//Read in simulation settings from .settings file
//Allows to change settings without recompiling
//constants.hpp still used for default settings
void PlasmaDomain::readSettingsFile(const char* settings_filename)
{
  std::ifstream in_file(settings_filename);
  std::string line;
  while(std::getline(in_file, line)){
    if(line[0] == '#') continue; //skip comment lines
    std::istringstream ss_line(line);
    std::string lhs, rhs;
    std::getline(ss_line,lhs,'=');
    std::getline(ss_line,rhs,';');
    auto it = std::find(m_setting_names.begin(),m_setting_names.end(),lhs);
    assert(it != m_setting_names.end());
    auto index = std::distance(m_setting_names.begin(),it);
    switch (static_cast<int>(index)) {
    case static_cast<int>(Setting::x_bound_1): x_bound_1 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Setting::x_bound_2): x_bound_2 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Setting::y_bound_1): y_bound_1 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Setting::y_bound_2): y_bound_2 = stringToBoundaryCondition(rhs); break;
    case static_cast<int>(Setting::radiative_losses): radiative_losses = (rhs == "true"); break;
    case static_cast<int>(Setting::ambient_heating): ambient_heating = (rhs == "true"); break;
    case static_cast<int>(Setting::thermal_conduction): thermal_conduction = (rhs == "true"); break;
    case static_cast<int>(Setting::flux_saturation): flux_saturation = (rhs == "true"); break;
    case static_cast<int>(Setting::temp_chromosphere): temp_chromosphere = std::stod(rhs); break;
    case static_cast<int>(Setting::radiation_ramp): radiation_ramp = std::stod(rhs); break;
    case static_cast<int>(Setting::heating_rate): heating_rate = std::stod(rhs); break;
    case static_cast<int>(Setting::b_0): b_0 = std::stod(rhs); break;
    case static_cast<int>(Setting::epsilon): epsilon = std::stod(rhs); break;
    case static_cast<int>(Setting::epsilon_thermal): epsilon_thermal = std::stod(rhs); break;
    case static_cast<int>(Setting::epsilon_rad): epsilon_rad = std::stod(rhs); break;
    case static_cast<int>(Setting::epsilon_viscous): epsilon_viscous = std::stod(rhs); break;
    case static_cast<int>(Setting::dt_thermal_min): dt_thermal_min = std::stod(rhs); break;
    case static_cast<int>(Setting::rho_min): rho_min = std::stod(rhs); break;
    case static_cast<int>(Setting::temp_min): temp_min = std::stod(rhs); break;
    case static_cast<int>(Setting::thermal_energy_min): thermal_energy_min = std::stod(rhs); break;
    case static_cast<int>(Setting::max_iterations): max_iterations = std::stoi(rhs); break;
    case static_cast<int>(Setting::max_time): max_time = std::stod(rhs); break;
    case static_cast<int>(Setting::iter_output_interval): iter_output_interval = std::stoi(rhs); break;
    case static_cast<int>(Setting::time_output_interval): time_output_interval = std::stod(rhs); break;
    case static_cast<int>(Setting::output_flags): /*TODO*/ break;
    case static_cast<int>(Setting::state_flags): /*TODO*/ break;
    case static_cast<int>(Setting::xdim): m_xdim = std::stoi(rhs); break;
    case static_cast<int>(Setting::ydim): m_ydim = std::stoi(rhs); break;
    case static_cast<int>(Setting::dx): m_dx = std::stod(rhs); break;
    case static_cast<int>(Setting::dy): m_dy = std::stod(rhs); break;
    default: break;
    }
  }
}

void PlasmaDomain::outputPreamble()
{
  m_out_file << m_xdim << "," << m_ydim << std::endl;
  m_out_file << m_dx << "," << m_dy << std::endl;
  m_out_file << m_grids[b_x].format(',',';');
  m_out_file << m_grids[b_y].format(',',';');
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
//All variables marked in m_state_flags are output here
void PlasmaDomain::writeStateFile(int precision) const
{
  std::ofstream state_file;
  state_file.open(m_run_name+std::to_string(m_iter%2)+".state");
  state_file << m_xdim << "," << m_ydim << std::endl;
  state_file << m_dx << "," << m_dy << std::endl;
  state_file << "t=" << m_time << std::endl;
  for(int i=0; i<num_variables; i++){
    if(m_state_flags[i]){
      state_file << m_var_names[i] << std::endl;
      state_file << m_grids[i].format(',',';',precision);
    }
  }
  state_file.close();
}

//Gets rid of the older of the two state files at the end of the sim run
void PlasmaDomain::cleanUpStateFiles() const
{
  rename((m_run_name+std::to_string((m_iter-1)%2)+".state").c_str(),
          (m_run_name+".state").c_str());
  remove((m_run_name+std::to_string((m_iter-2)%2)+".state").c_str());
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
  printf("Iter: %i",m_iter);
  if(max_iterations > 0) printf("/%i",max_iterations);
  printf("|t: %f",m_time);
  if(max_time > 0.0) printf("/%f",max_time);
  printf("|dt: %f",min_dt);
  if(thermal_conduction) printf("|Cond. Subcyc: %i",subcycles_thermal);
  if(radiative_losses) printf("|Rad. Subcyc: %i",subcycles_rad);
  printf("\n");
}