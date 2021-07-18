//fileio.cpp
//PlasmaDomain functionality relating to file I/O

#include <string>
#include <fstream>
#include <sstream>
#include <cstdio>
#include "plasmadomain.hpp"

//Corresponding variable names for file I/O
//Must match the ordering of the Variable enum defined in plasmadomain.hpp
const std::vector<std::string> PlasmaDomain::m_var_names = {
  "rho","temp","mom_x","mom_y","b_x","b_y","b_z",
  "press","energy","rad","dt","dt_thermal","dt_rad","v_x","v_y",
  "grav_x","grav_y","b_magnitude","b_hat_x","b_hat_y",
  "mag_press","mag_pxx","mag_pyy","mag_pzz","mag_pxy","mag_pxz","mag_pyz"
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
void PlasmaDomain::writeStateFile() const
{
  std::ofstream state_file;
  state_file.open(m_run_name+std::to_string(m_iter%2)+".state");
  state_file << m_xdim << "," << m_ydim << std::endl;
  state_file << m_dx << "," << m_dy << std::endl;
  state_file << "t=" << m_time << std::endl;
  for(int i=0; i<num_variables; i++){
    if(m_state_flags[i]){
      state_file << m_var_names[i] << std::endl;
      state_file << m_grids[i].format(',',';',-1);
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
