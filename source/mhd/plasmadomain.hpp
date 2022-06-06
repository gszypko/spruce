//plasmadomain.hpp
//Defines PlasmaDomain class for containing parameters
//and time evolution of a 2D domain of plasma

#ifndef PLASMADOMAIN_HPP
#define PLASMADOMAIN_HPP

#include <filesystem>
namespace fs = std::filesystem;
#include <vector>
#include <string>
#include "grid.hpp"
#include "modulehandler.hpp"
#include "equationset.hpp"

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

  //Constructors and Initialization
  PlasmaDomain(const fs::path &out_path, const fs::path &config_path, const fs::path &state_file, bool continue_mode, bool overwrite_init);
  PlasmaDomain() = default;
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
  friend class FieldHeating;
  friend class RadiativeLosses;
  friend class ThermalConduction;

  //Friend class declarations (for EquationSets)
  friend class EquationSet;
  friend class IdealMHD;

  enum InternalVars {d_x,d_y,pos_x,pos_y,be_x,be_y};
  static const inline std::vector<std::string> m_internal_var_names = {"d_x","d_y","pos_x","pos_y","be_x","be_y"};
  std::vector<Grid> m_internal_grids{m_internal_var_names.size()};
  Grid &m_d_x = m_internal_grids[d_x], &m_d_y = m_internal_grids[d_y],
  &m_pos_x = m_internal_grids[pos_x], &m_pos_y  = m_internal_grids[pos_y],
  &m_be_x = m_internal_grids[be_x], &m_be_y = m_internal_grids[be_y];

  //Strings corresponding to variables, settings, boundary conditions for file I/O
  int m_xl, m_xu, m_yl, m_yu; //Lower and upper bounds for diff'l operations on the domain (excluding ghost zones)
  Grid m_ghost_zone_mask; //Equals 0 inside ghost zones and 1 everywhere else; for multiplying to negate values in ghost zones

  Grid& grid(int var) { return m_eqs->grid(var); }
  // std::vector<bool> m_output_flags; //Variables that are printed in .out files (for visualization purposes)
  double m_ion_mass; //in g
  double m_adiabatic_index; //aka gamma, unitless

  ModuleHandler m_module_handler;
  std::unique_ptr<EquationSet> m_eqs;

  fs::path m_out_directory;
  std::ofstream m_out_file;
  bool m_overwrite_init;
  std::vector<std::string> comment_lines;

  double m_time;
  double m_duration;
  int m_iter;
  double max_time; //Upper bound on simulation time

  /**************************** CONFIGS ******************************/
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
 
  enum class Config {
    x_bound_1, x_bound_2, y_bound_1, y_bound_2,
    epsilon, epsilon_viscous, rho_min,
    temp_min, thermal_energy_min, max_iterations, iter_output_interval, time_output_interval,
    output_flags, xdim, ydim, open_boundary_strength, std_out_interval, safe_state_mode, safe_state_interval,
    open_boundary_decay_base, x_origin, y_origin, time_integrator, equation_set
  };
  static inline std::vector<std::string> m_config_names = {
    "x_bound_1","x_bound_2","y_bound_1","y_bound_2",
    "epsilon","epsilon_viscous","rho_min",
    "temp_min","thermal_energy_min","max_iterations","iter_output_interval","time_output_interval",
    "output_flags","xdim","ydim","open_boundary_strength","std_out_interval","safe_state_mode", "safe_state_interval",
    "open_boundary_decay_base", "x_origin", "y_origin", "time_integrator", "equation_set"
  };

  /*********************************************************************/ 


  void advanceTime(bool verbose = true);

  void integrateEuler(double time_step, double visc_coeff);
  void integrateRK2(double time_step, double visc_coeff);
  void integrateRK4(double time_step, double visc_coeff);

  void updateGhostZones();

  void openBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);
  void reflectBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);
  void fixedBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);

  // void recomputeDT();

  // void computeConstantTerms();
  void computeIterationBounds();

  void outputPreamble();
  void outputCurrentState();
  void writeGridToOutput(const Grid& grid, std::string var_name);
  void writeStateFile(std::string filename_stem = "mhd",int precision = -1) const;
  void cleanUpStateFiles(std::string filename_stem = "end") const;
  void printUpdate(double dt) const;

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
  
  double boundaryInterpolate(const Grid &quantity, int i1, int j1, int i2, int j2);
  double boundaryExtrapolate(const Grid &quantity, int i1, int j1, int i2, int j2);

};

#endif