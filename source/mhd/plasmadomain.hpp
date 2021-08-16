//plasmadomain.hpp
//Defines PlasmaDomain class for containing parameters
//and time evolution of a 2D domain of plasma

#ifndef PLASMADOMAIN_HPP
#define PLASMADOMAIN_HPP

#include "utils.hpp"
#include "grid.hpp"
#include "constants.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <filesystem>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <omp.h>

class PlasmaDomain
{
public:
  enum class BoundaryCondition { Periodic, Wall, Open };

  //state variables (taken as input, carried over between iterations)
  enum StateVars { state_var_start=0,
                        pos_x=state_var_start,pos_y,rho,temp,mom_x,mom_y,b_x,b_y,b_z,grav_x,grav_y,
                        state_var_end };
  
  //derived variables (derived from state variables)
  enum DerivedVars { derived_var_start=state_var_end,
                          press=derived_var_start,energy,rad,dt,dt_thermal,dt_rad,v_x,v_y,
                          derived_var_end };

  //constant variables (unchanging bw iterations, derived from state variables)
  enum ConstantVars { constant_var_start=derived_var_end,
                           b_magnitude=constant_var_start,b_hat_x,b_hat_y,d_x,d_y,mag_press,mag_pxx,mag_pyy,mag_pzz,mag_pxy,mag_pxz,mag_pyz,
                           constant_var_end, num_variables=constant_var_end };

  //Constructors and Initialization
  PlasmaDomain(const char* out_path, const char* config_path, bool continue_mode = false);
  void initialize(const std::vector<Grid>& input_vars, double ion_mass, double adiabatic_index);
  void hydrostaticInitialize();
  void gaussianInitialize(double min_rho, double max_rho, double min_temp, double max_temp, double std_dev_x, double std_dev_y);

  void readStateFile(const char* in_filename);
  void readConfigFile(const char* settings_filename);

  //Time Evolution
  void run(double time_duration);

private:
  //Strings corresponding to variables, settings, boundary conditions for file I/O
  //TODO: Move initialization here, now that no longer static!
  static const std::vector<std::string> m_var_names;
  static const std::vector<std::string> m_config_names;
  static const std::vector<std::string> m_boundary_condition_names;
  int m_xl, m_xu, m_yl, m_yu; //Lower and upper bounds for diff'l operations on the domain (excluding ghost zones)

  std::vector<Grid> m_grids;
  std::vector<bool> m_output_flags; //Variables that are printed in .out files (for visualization purposes)
  double m_ion_mass; //in g
  double m_adiabatic_index; //aka gamma, unitless

  std::filesystem::path m_out_directory;
  std::ofstream m_out_file;

  double m_time;
  int m_iter;
  double max_time; //Upper bound on simulation time

  /**************************** SETTINGS ******************************/
  size_t m_xdim, m_ydim;
  int max_iterations; //Upper bound on simulation iterations; unbounded if negative
  //Boundary condition settings
  BoundaryCondition x_bound_1, x_bound_2, y_bound_1, y_bound_2;
  //Physics settings
  bool radiative_losses, ambient_heating, thermal_conduction;
  bool flux_saturation; //Config for thermal conduction
  //Physical parameters
  double temp_chromosphere; //Cutoff temperature for radiative losses
  double radiation_ramp;  //Width of cutoff ramp, in units of temperature, for low-temp radiation
  double heating_rate;  //Constant ambient heating rate
  //Safety factors
  double epsilon, epsilon_thermal, epsilon_rad; //Time step calculation
  double epsilon_viscous; //Prefactor for artificial viscosity
  double dt_thermal_min; //Minimum timestep for thermal conduction
  double rho_min, temp_min, thermal_energy_min; //Lower bounds for mass density and thermal energy density
  //Output settings
  int iter_output_interval;
  double time_output_interval;
  bool continue_mode; //true when continuing previous run; appends results to mhd.out and replaces mhd.state
  /*********************************************************************/ 

  void advanceTime(bool verbose = true);
  void subcycleConduction(int subcycles_thermal, double dt_total);
  void subcycleRadiation(int subcycles_rad, double dt_total);

  void recomputeDerivedVariables();
  void recomputeTemperature();
  void computeConstantTerms();
  void computeIterationBounds();

  void catchUnderdensity();
  void updateGhostZones();
  void openBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);
  void wallBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);

  void recomputeRadiativeLosses();
  void recomputeDT();
  void recomputeDTThermal();
  void recomputeDTRadiative();
  
  void outputPreamble();
  void outputCurrentState();
  void writeStateFile(int precision = -1) const;
  void cleanUpStateFiles() const;
  void printUpdate(double min_dt, int subcycles_thermal, int subcycles_rad) const;

  void setOutputFlag(std::string var_name, bool new_flag);
  void setOutputFlags(const std::vector<std::string> var_names, bool new_flag);

  BoundaryCondition stringToBoundaryCondition(const std::string str) const;
  void handleSingleConfig(int setting_index, std::string rhs);
  void handleConfigList(int setting_index, std::vector<std::string> rhs_vec,
                        std::vector<std::vector<std::string> > &rhs_lists,
                        std::vector<int> &list_vars, int &num_combinations);

  //Compute surface values from cell-centered values using Barton's method
  //Meant to be used for transport terms only
  //Result indexed s.t. element i,j indicates surface between i,j and i-1,j
  //if "index"==0, or i,j and i,j-1 if "index"==1
  Grid upwindSurface(const Grid &cell_center, const Grid &vel, const int index);

  Grid transportDerivative1D(const Grid &quantity, const Grid &vel, const int index);

  //Compute single-direction divergence term for non-transport term (central differencing)
  Grid derivative1D(const Grid &quantity, const int index);

  //Compute divergence term for simulation parameter "quantity"
  //"quantity","vx","vy" used for transport term
  //Non-transport terms contained in "nontransp_x", "nontransp_y"
  Grid divergence(const Grid &quantity, const Grid &nontransp_x, const Grid &nontransp_y);

  //Compute single-direction second derivative
  Grid secondDerivative1D(const Grid &quantity, const int index);

  //Computes Laplacian (del squared) of "quantity"
  Grid laplacian(const Grid &quantity);

  //Computes 1D cell-centered conductive flux from temperature "temp"
  //Flux computed in direction indicated by "index": 0 for x, 1 for y
  //k0 is conductive coefficient
  Grid oneDimConductiveFlux(const Grid &temp, const Grid &rho, double k0, int index);

  //Computes cell-centered, field-aligned conductive flux from temperature "temp"
  //temp is temperature Grid
  //b_hat_x, b_hat_y are the components of the *unit* vector b_hat
  //k0 is conductive coefficient
  //Output is written to flux_out_x and flux_out_y
  void fieldAlignedConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &temp, const Grid &rho,
                                      const Grid &b_hat_x, const Grid &b_hat_y, const double k0);

  //Computes saturated conductive flux at each point in grid,
  //then ensures that provided fluxes do not exceed the saturation point
  void saturateConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &rho, const Grid &temp);
};

#endif