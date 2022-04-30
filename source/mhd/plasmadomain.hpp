//plasmadomain.hpp
//Defines PlasmaDomain class for containing parameters
//and time evolution of a 2D domain of plasma

#ifndef PLASMADOMAIN_HPP
#define PLASMADOMAIN_HPP

#include "utils.hpp"
#include "grid.hpp"
#include "constants.hpp"
#include "modulehandler.hpp"
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <filesystem>
namespace fs = std::filesystem;
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <omp.h>

class PlasmaDomain
{
public:
  enum class BoundaryCondition { Periodic, Open, Fixed, Reflect };
  static inline std::vector<std::string> m_boundary_condition_names = {
    "periodic", "open", "fixed", "reflect"
  };
  enum class TimeIntegrator { Euler, RK2, RK4 };
  static inline std::vector<std::string> m_time_integrator_names = {
    "euler", "rk2", "rk4"
  };

  enum Vars { d_x,d_y,pos_x,pos_y,rho,temp,mom_x,mom_y,be_x,be_y,grav_x,grav_y,
              press,thermal_energy,kinetic_energy,div_bi,bi_x,bi_y,
              dt,v_x,v_y,n,
              div_be,b_hat_x,b_hat_y,b_magnitude,v_a,
              num_variables };
  static inline std::vector<std::string> m_var_names = {
    "d_x","d_y","pos_x","pos_y","rho","temp","mom_x","mom_y","be_x","be_y","grav_x","grav_y",
    "press","thermal_energy","kinetic_energy","div_bi","bi_x","bi_y","dt","v_x","v_y","n",
    "div_be","b_hat_x","b_hat_y","b_magnitude","v_a"
  };
  static inline std::vector<int> state_vars = {d_x,d_y,pos_x,pos_y,rho,temp,mom_x,mom_y,be_x,be_y,bi_x,bi_y,grav_x,grav_y};

  //Constructors and Initialization
  PlasmaDomain(const fs::path &out_path, const fs::path &config_path, const std::vector<Grid>& input_vars, double ion_mass, double adiabatic_index);
  PlasmaDomain(const fs::path &out_path, const fs::path &config_path, const fs::path &state_file, bool continue_mode);
  void readStateFile(const fs::path &state_file, bool continue_mode = true);
  void readConfigFile(const fs::path &config_file);

  //Time Evolution
  void run(double time_duration);

  //Takes not-necessarily-uniform cell sizes d_x and d_y and converts to cell center positions,
  //returned as the first (x position) and second (y position) indicies of a vector of Grids
  //x_origin and y_origin specify location of the origin (0,0) in the domain:
  //"lower" is the default, "center" and "upper" are also options for each
  //For example, x_origin = "center" and y_origin = "center" places the origin
  //at the center of the domain.
  static Grid convertCellSizesToCellPositions(const Grid& d, int index, std::string origin_position);

  static bool validateCellSizesAndPositions(const Grid& d, const Grid& pos, int index, double tolerance = 1.0e-4);

private:

  //Friend class declarations (for Modules)
  friend class ModuleHandler;
  friend class Module;
  friend class SGFilter;
  friend class AmbientHeating;
  friend class LocalizedHeating;
  friend class RadiativeLosses;
  friend class ThermalConduction;

  //Strings corresponding to variables, settings, boundary conditions for file I/O
  int m_xl, m_xu, m_yl, m_yu; //Lower and upper bounds for diff'l operations on the domain (excluding ghost zones)
  Grid m_ghost_zone_mask; //Equals 0 inside ghost zones and 1 everywhere else; for multiplying to negate values in ghost zones

  std::vector<Grid> m_grids;
  std::vector<bool> m_output_flags; //Variables that are printed in .out files (for visualization purposes)
  double m_ion_mass; //in g
  double m_adiabatic_index; //aka gamma, unitless

  ModuleHandler m_module_handler;

  fs::path m_out_directory;
  std::ofstream m_out_file;

  double m_time;
  int m_iter;
  double max_time; //Upper bound on simulation time

  /**************************** SETTINGS ******************************/
  size_t m_xdim, m_ydim;
  int max_iterations; //Upper bound on simulation iterations; unbounded if negative
  //Boundary condition settings
  BoundaryCondition x_bound_1, x_bound_2, y_bound_1, y_bound_2;
  double open_boundary_strength; // multiple of local sound speed added to velocity at boundary surface
  double open_boundary_decay_base; // base of exponential decay of rho, thermal_energy beyond open boundary surface
  //Safety factors
  double epsilon; //Time step calculation
  double epsilon_viscous; //Prefactor for artificial viscosity
  double rho_min, temp_min, thermal_energy_min; //Lower bounds for mass density and thermal energy density
  TimeIntegrator time_integrator; //indicates time integration scheme to use
  //Output settings
  int iter_output_interval;
  double time_output_interval;
  int std_out_interval; //number of iterations between printing an update to standard out
  bool safe_state_mode; //when true, outputs a state file for every iteration; when false, only outputs a state when the run completes without issue
  int safe_state_interval; //when safe_state_mode == true, number of iterations between .state files written out during run
  bool continue_mode; //true when continuing previous run; appends results to mhd.out and replaces mhd.state
  // std::string x_origin, y_origin; //specifies where (0,0) position is located; each can be one of "lower", "center", or "upper"
 
  enum class Config {
    x_bound_1, x_bound_2, y_bound_1, y_bound_2,
    epsilon, epsilon_viscous, rho_min,
    temp_min, thermal_energy_min, max_iterations, iter_output_interval, time_output_interval,
    output_flags, xdim, ydim, open_boundary_strength, std_out_interval, safe_state_mode, safe_state_interval,
    open_boundary_decay_base, x_origin, y_origin, time_integrator
  };
  static inline std::vector<std::string> m_config_names = {
    "x_bound_1","x_bound_2","y_bound_1","y_bound_2",
    "epsilon","epsilon_viscous","rho_min",
    "temp_min","thermal_energy_min","max_iterations","iter_output_interval","time_output_interval",
    "output_flags","xdim","ydim","open_boundary_strength","std_out_interval","safe_state_mode", "safe_state_interval",
    "open_boundary_decay_base", "x_origin", "y_origin", "time_integrator"
  };

  /*********************************************************************/ 

  void configureSimulation(const fs::path &config_path, const fs::path &out_path, bool continue_mode = false);

  void advanceTime(bool verbose = true);

  void integrateEuler(std::vector<Grid> &grids, double time_step, double visc_coeff);
  void integrateRK2(std::vector<Grid> &grids, double time_step, double visc_coeff);
  void integrateRK4(std::vector<Grid> &grids, double time_step, double visc_coeff);
  std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff);
  void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step);

  //Functions to be applied between time steps
  void propagateChanges(); //Encapsulates all of the below, except SG filtering
  void recomputeDerivedVariables();
  void recomputeTemperature();
  void catchUnderdensity();
  void updateGhostZones();

  //Versions of the above for operating on intermediate states of the system
  //i.e. for acting on Grid vectors other than m_grids
  void propagateChanges(std::vector<Grid> &grids);
  void recomputeDerivedVariables(std::vector<Grid> &grids);
  void recomputeTemperature(std::vector<Grid> &grids);
  void catchUnderdensity(std::vector<Grid> &grids);
  void updateGhostZones(std::vector<Grid> &grids);

  void openBoundaryExtrapolate(std::vector<Grid> &grids, int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);
  void reflectBoundaryExtrapolate(std::vector<Grid> &grids, int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);
  void fixedBoundaryExtrapolate(std::vector<Grid> &grids, int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);

  void recomputeDT();

  void computeConstantTerms();
  void computeIterationBounds();

  void outputPreamble();
  void outputCurrentState();
  void writeGridToOutput(const Grid& grid, std::string var_name);
  void writeStateFile(std::string filename_stem = "mhd",int precision = -1) const;
  void cleanUpStateFiles(std::string filename_stem = "end") const;
  void printUpdate(double dt) const;

  void setOutputFlag(std::string var_name, bool new_flag);
  void setOutputFlags(const std::vector<std::string> var_names, bool new_flag);

  BoundaryCondition stringToBoundaryCondition(const std::string str) const;
  TimeIntegrator stringToTimeIntegrator(const std::string str) const;
  void handleSingleConfig(int setting_index, std::string rhs);
  void handleConfigList(int setting_index, std::vector<std::string> rhs_vec,
                        std::vector<std::vector<std::string> > &rhs_lists,
                        std::vector<int> &list_vars, int &num_combinations);

  //Compute surface values from cell-centered values using Barton's method
  //Meant to be used for transport terms only
  //Result indexed s.t. element i,j indicates surface between i,j and i-1,j
  //if "index"==0, or i,j and i,j-1 if "index"==1
  Grid upwindSurface(const Grid &cell_center, const Grid &vel, const int index);

  //Compute divergence of scalar quantity times velocity using Barton's method
  //Meant for advection terms
  Grid transportDivergence2D(const Grid &quantity, const std::vector<Grid> &vel);
  //Compute 1D derivative using Barton's method, for transport terms
  Grid transportDerivative1D(const Grid &quantity, const Grid &vel, const int index);

  //Compute divergence of vector quantity a_x,a_y
  Grid divergence2D(const Grid& a_x, const Grid& a_y);
  //Compute single-direction divergence term for non-transport term (central differencing)
  Grid derivative1D(const Grid &quantity, const int index);

  //Compute divergence of vector quantity a[0],a[1]
  Grid divergence2D(const std::vector<Grid>& a);

  //Compute single-direction second derivative
  Grid secondDerivative1D(const Grid &quantity, const int index);

  //Computes Laplacian (del squared) of "quantity"
  Grid laplacian(const Grid &quantity);

  //Computes curl of vector in z-direction (result in xy-plane)
  std::vector<Grid> curlZ(const Grid& z);

  //Computes curl of vector in xy-plane (result in z-direction)
  Grid curl2D(const Grid& x, const Grid& y);
  
  //Compute (perturbation) magnetic field term on RHS of energy equation
  Grid computeMagneticEnergyTerm();

  double boundaryInterpolate(const Grid &quantity, int i1, int j1, int i2, int j2);
  double boundaryExtrapolate(const Grid &quantity, int i1, int j1, int i2, int j2);

};

#endif