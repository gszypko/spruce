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
  enum BoundaryCondition { Periodic, Wall, Open };

  /******** SETTINGS (to be accessed and set directly, currently) ********/
  /*** All settings false or zero by default unless modified after construction ***/
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
  int chromosphere_depth; //In number of Grid cells
  //Safety factors
  double epsilon, epsilon_thermal, epsilon_rad; //Time step calculation
  double epsilon_viscous; //Prefactor for artificial viscosity
  double dt_thermal_min; //Minimum timestep for thermal conduction
  double rho_min, thermal_energy_min; //Lower bounds for mass density and thermal energy density
  //Output settings
  int n_iterations, output_interval;
  /*********************************************************************/ 

  PlasmaDomain(size_t xdim, size_t ydim, double dx, double dy, const char* run_name);
  PlasmaDomain(const char* run_name);
  PlasmaDomain();

  //Implemented
  void readStateFile(const char* in_filename);
  void setDefaultSettings();
  void hydrostaticInitialize();
  void gaussianInitialize();
  void outputPreamble();
  void outputCurrentState();
  void setSolarGravity(double base_gravity, double r_solar);
  void writeStateFile() const;
  void setOutputFlag(int var, bool new_flag = true);
  void setOutputFlags(const std::vector<int> vars, bool new_flag = true);
  void setStateFlag(int var, bool new_flag = true);
  void setStateFlags(const std::vector<int> vars, bool new_flag = true);
  void advanceTime(bool verbose = true);
  void cleanUpStateFiles() const;

  //Not yet implemented
  // void clampWallBoundaries();
  // void cleanUpStateFiles();

private:
  void recomputeDerivedVariables();
  void computeMagneticTerms();
  void recomputeTemperature();
  void catchUnderdensity();
  void clampWallBoundaries(Grid& mom_x_next, Grid& mom_y_next, Grid& rho_next, Grid& energy_next);
  void recomputeRadiativeLosses();
  void recomputeDT();
  void recomputeDTThermal();
  void recomputeDTRadiative();
  void subcycleConduction(int subcycles_thermal, double dt_total);
  void subcycleRadiation(int subcycles_rad, double dt_total);

  static const std::vector<std::string> m_var_names;
  std::vector<Grid> m_grids;
  std::vector<bool> m_output_flags; //Variables that are printed in .out files (for visualization purposes)
  std::vector<bool> m_state_flags; //Variables that are printed in .state files (should be a minimal complete description of the plasma)

  size_t m_xdim, m_ydim;
  double m_dx, m_dy;

  std::string m_run_name;
  std::ofstream m_out_file;

  double m_time;
  int m_iter;
};

#endif