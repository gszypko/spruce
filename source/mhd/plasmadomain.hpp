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
#include <filesystem>
#include "grid.hpp"

class PlasmaDomain
{
public:
  enum class BoundaryCondition { Periodic, Wall, Open };
  enum Variable {
    pos_x,pos_y,rho,temp,mom_x,mom_y,b_x,b_y,b_z,grav_x,grav_y, //base variables (taken as input, carried over between iterations)
    press,energy,rad,dt,dt_thermal,dt_rad,v_x,v_y, //derived variables (derived from base variables)
    b_magnitude,b_hat_x,b_hat_y,d_x,d_y, //constant variables (unchanging bw iterations)
    mag_press,mag_pxx,mag_pyy,mag_pzz,mag_pxy,mag_pxz,mag_pyz, //constant variables
    num_variables //never add new variable after this in the enum!
  };

  //Constructors and Initialization
  PlasmaDomain(const char* run_name, const char* settings_file_name);
  void initialize(const Grid& rho, const Grid& temp, const Grid& mom_x, const Grid& mom_y,
                  const Grid& b_x, const Grid& b_y, const Grid& b_z,
                  const Grid& pos_x, const Grid& pos_y, const Grid& grav_x, const Grid& grav_y);
  void hydrostaticInitialize();
  void gaussianInitialize(double min_rho, double max_rho, double min_temp, double max_temp, double std_dev_x, double std_dev_y);

  void readStateFile(const char* in_filename);
  void readConfigFile(const char* settings_filename);

  //Time Evolution
  void run();

private:
  //Strings corresponding to variables, settings, boundary conditions for file I/O
  static const std::vector<std::string> m_var_names;
  static const std::vector<std::string> m_config_names;
  static const std::vector<std::string> m_boundary_condition_names;
  int m_xl, m_xu, m_yl, m_yu; //Lower and upper bounds for diff'l operations on the domain (excluding ghost zones)

  std::vector<Grid> m_grids;
  std::vector<bool> m_output_flags; //Variables that are printed in .out files (for visualization purposes)
  std::vector<bool> m_state_flags; //Variables that are printed in .state files (should be a minimal complete description of the plasma)

  std::filesystem::path m_out_directory;
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
  bool flux_saturation; //Config for thermal conduction
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
  void setStateFlag(std::string var_name, bool new_flag);
  void setStateFlags(const std::vector<std::string> var_names, bool new_flag);

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