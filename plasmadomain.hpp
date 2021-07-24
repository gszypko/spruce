//plasmadomain.hpp
//Defines PlasmaDomain class for containing parameters
//and time evolution of a 2D domain of plasma

#ifndef PLASMADOMAIN_HPP
#define PLASMADOMAIN_HPP

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <map>
#include <iostream>
#include "grid.hpp"

class PlasmaDomain
{
public:
  enum class BoundaryCondition { Periodic, Wall, Open };
  enum Variable {
    rho,temp,mom_x,mom_y,b_x,b_y,b_z, //base variables (carried over between iterations)
    press,energy,rad,dt,dt_thermal,dt_rad,v_x,v_y, //derived variables (derived from base variables)
    grav_x,grav_y,b_magnitude,b_hat_x,b_hat_y, //constant variables (unchanging bw iterations)
    mag_press,mag_pxx,mag_pyy,mag_pzz,mag_pxy,mag_pxz,mag_pyz, //constant variables
    num_variables //never add new variable after this in the enum!
  };

  //Constructors and Initialization
  PlasmaDomain(size_t xdim, size_t ydim, double dx, double dy, const char* run_name);
  PlasmaDomain(const char* run_name, const char* settings_file_name);
  PlasmaDomain(const char* run_name);
  PlasmaDomain();
  void setDefaultSettings();
  void hydrostaticInitialize();
  void gaussianInitialize(double min_rho, double max_rho, double min_temp, double max_temp, double std_dev_x, double std_dev_y);
  void setSolarGravity(double base_gravity, double r_solar);

  void readStateFile(const char* in_filename);
  void readSettingsFile(const char* settings_filename);

  //Time Evolution
  void run();

private:
  //Strings corresponding to variables, settings, boundary conditions for file I/O
  static const std::vector<std::string> m_var_names;
  static const std::vector<std::string> m_setting_names;
  static const std::vector<std::string> m_boundary_condition_names;

  std::vector<Grid> m_grids;
  std::vector<bool> m_output_flags; //Variables that are printed in .out files (for visualization purposes)
  std::vector<bool> m_state_flags; //Variables that are printed in .state files (should be a minimal complete description of the plasma)

  std::string m_run_name;
  std::ofstream m_out_file;

  double m_time;
  int m_iter;

  /**************************** SETTINGS ******************************/
  size_t m_xdim, m_ydim;
  double m_dx, m_dy;
  int max_iterations; //Upper bound on simulation iterations; unbounded if negative
  double max_time; //Upper bound on simulation time; unbounded if negative
  //Boundary condition settings
  BoundaryCondition x_bound_1, x_bound_2, y_bound_1, y_bound_2;
  //Physics settings
  bool radiative_losses, ambient_heating, thermal_conduction;
  bool flux_saturation; //Setting for thermal conduction
  //Physical parameters
  double temp_chromosphere; //Minimum allowed temperature
  double radiation_ramp;  //Width of cutoff ramp, in units of temperature, for low-temp radiation
  double heating_rate;  //Constant ambient heating rate
  double b_0;  //Base value of magnetic field
  //Safety factors
  double epsilon, epsilon_thermal, epsilon_rad; //Time step calculation
  double epsilon_viscous; //Prefactor for artificial viscosity
  double dt_thermal_min; //Minimum timestep for thermal conduction
  double rho_min, temp_min, thermal_energy_min; //Lower bounds for mass density and thermal energy density
  //Output settings
  int iter_output_interval;
  double time_output_interval;
  /*********************************************************************/ 

  void advanceTime(bool verbose = true);
  void subcycleConduction(int subcycles_thermal, double dt_total);
  void subcycleRadiation(int subcycles_rad, double dt_total);

  void recomputeDerivedVariables();
  void recomputeTemperature();
  void computeMagneticTerms();

  void catchUnderdensity();
  void clampWallBoundaries(Grid& mom_x_next, Grid& mom_y_next, Grid& rho_next, Grid& energy_next);

  void recomputeRadiativeLosses();
  void recomputeDT();
  void recomputeDTThermal();
  void recomputeDTRadiative();

  void outputPreamble();
  void outputCurrentState();
  void writeStateFile(int precision = -1) const;
  void cleanUpStateFiles() const;
  void printUpdate(double min_dt, int subcycles_thermal, int subcycles_rad) const;

  void setOutputFlag(int var, bool new_flag = true);
  void setOutputFlags(const std::vector<int> vars, bool new_flag = true);
  void setStateFlag(int var, bool new_flag = true);
  void setStateFlags(const std::vector<int> vars, bool new_flag = true);

  BoundaryCondition stringToBoundaryCondition(const std::string str) const;
};

#endif