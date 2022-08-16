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
#include "savitzkygolay.hpp"

class PlasmaDomain
{
public: //******************************************************************************************************
  // Boundary Condition Options
  enum class BoundaryCondition { Periodic, Open, Fixed, Reflect, OpenMoC };
  static inline std::vector<std::string> m_boundary_condition_names = {
    "periodic", "open", "fixed", "reflect", "open_moc"
  };

  // Time Integrator Options
  enum class TimeIntegrator { Euler, RK2, RK4 };
  static inline std::vector<std::string> m_time_integrator_names = {
    "euler", "rk2", "rk4"
  };

  // Internal Grids
  enum Grids {d_x,d_y,pos_x,pos_y,be_x,be_y};
  static const inline std::vector<std::string> m_gridnames = {"d_x","d_y","pos_x","pos_y","be_x","be_y"};
  std::vector<Grid> m_grids{m_gridnames.size()};

  // Config Names
  enum class Config {
    x_bound_1, x_bound_2, y_bound_1, y_bound_2,
    epsilon, epsilon_viscous, density_min,
    temp_min, thermal_energy_min, max_iterations, iter_output_interval, time_output_interval,
    output_flags, xdim, ydim, open_boundary_strength, std_out_interval, write_interval,
    open_boundary_decay_base, x_origin, y_origin, time_integrator, equation_set, duration, epsilon_courant
  };
  static inline std::vector<std::string> m_config_names = {
    "x_bound_1","x_bound_2","y_bound_1","y_bound_2",
    "epsilon","epsilon_viscous","density_min",
    "temp_min","thermal_energy_min","max_iterations","iter_output_interval","time_output_interval",
    "output_flags","xdim","ydim","open_boundary_strength","std_out_interval","write_interval",
    "open_boundary_decay_base", "x_origin", "y_origin", "time_integrator", "equation_set", "duration", "epsilon_courant"
  };

  // Public Member Functions
  PlasmaDomain(const fs::path &out_path, const fs::path &config_path, const fs::path &state_file, bool continue_mode, bool overwrite_init);
  PlasmaDomain(): m_module_handler{*this} {};
  void readStateFile(const fs::path &state_file, bool continue_mode = true);
  void readConfigFile(const fs::path &config_file);
  void initOutputContainer();
  static Grid convertCellSizesToCellPositions(const Grid& d, int index, std::string origin_position);
  static bool validateCellSizesAndPositions(const Grid& d, const Grid& pos, int index, double tolerance = 1.0e-4);
  void run(double time_duration);

private: //*****************************************************************************************************
  //Friend class declarations (for Modules)
  friend class ModuleHandler;
  friend class Module;
  friend class SGFilter;
  friend class AmbientHeating;
  friend class LocalizedHeating;
  friend class FieldHeating;
  friend class RadiativeLosses;
  friend class ThermalConduction;
  friend class CoulombExplosion;
  friend class EICThermalization;

  //Friend class declarations (for EquationSets)
  friend class EquationSet;
  friend class IdealMHD;
  friend class IdealMHDCons;
  friend class IdealMHD2E;
  friend class Ideal2F;

  // Time-Related Quantities
  TimeIntegrator m_time_integrator{TimeIntegrator::Euler}; //indicates time integration scheme to use
  double m_time{0.0};
  double m_duration{-1.0}; // 
  int m_iter{0}; // number of times the advanceTime() function has been called during simulation
  double m_max_time{-1.0}; // upper bound on simulation time
  int max_iterations{100}; // upper bound on simulation iterations; unbounded if negative

  // File Output Settings
  bool m_overwrite_init; // true when init.state to be overwritten at start of new simulation
  bool m_continue_mode; // true when continuing previous run; appends results to mhd.out and replaces mhd.state
  fs::path m_out_directory; // relative file path to where simulation files are to be saved
  fs::path m_out_filename{"mhd.out"}; // file name for saving time evolution of plasma
  std::vector<std::string> m_comment_lines{};
  int m_iter_output_interval{1};
  double m_time_output_interval{-1.0};
  int m_std_out_interval{1}; //number of iterations between printing an update to standard out
  int m_write_interval{1};
  int m_output_counter{0};
  std::vector<std::string> m_data_to_write; // holds data to be written to mhd.out
  int m_lines_recorded{0}; // holds the last line of m_data_to_write that was updated

  // variables relevant for boundary conditions
  int m_xl, m_xu, m_yl, m_yu; //Lower and upper bounds for diff'l operations on the domain (excluding ghost zones)
  Grid m_ghost_zone_mask; //Equals 0 inside ghost zones and 1 everywhere else; for multiplying to negate values in ghost zones
  BoundaryCondition x_bound_1, x_bound_2, y_bound_1, y_bound_2;
  double open_boundary_strength{0.0}; // multiple of local sound speed added to velocity at boundary surface
  double open_boundary_decay_base{1.0}; // base of exponential decay of rho, thermal_energy beyond open boundary surface

  // PlasmaDomain internal grids
  size_t m_xdim, m_ydim;
  Grid& grid(int index);
  Grid& grid(const std::string& name);
  int gridname2index(const std::string& name) const;
  bool is_grid(const std::string& name) const;
  
  // ion properties
  double m_ion_mass; //in g
  double m_adiabatic_index; //aka gamma, unitless

  // various class containers
  ModuleHandler m_module_handler;
  std::unique_ptr<EquationSet> m_eqs;
  SavitzkyGolay m_sg{};  
  
  // safety factors
  double epsilon; //Time step calculation
  double epsilon_viscous; //Prefactor for artificial viscosity
  double epsilon_courant; // Courant safety factor for 2F model
  double density_min, temp_min, thermal_energy_min; //Lower bounds for mass density and thermal energy density

  void advanceTime(bool verbose = true);

  void integrateEuler(double time_step, double visc_coeff);
  void integrateRK2(double time_step, double visc_coeff);
  void integrateRK4(double time_step, double visc_coeff);

  //******************** PRIVATE FUNCTION DECLARATIONS *************************************************************

  void updateGhostZones();

  void openBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);
  void reflectBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);
  void fixedBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4);

  void computeIterationBounds();
  bool allInternalGridsInitialized();

  void outputPreamble();
  void storeGrids();
  void writeToOutFile();
  void writeStateFile(std::string filename_stem = "mhd",int precision = -1) const;
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
  Grid upwindSurface(const Grid &cell_center, const Grid &vel, const int index) { return upwindSurface(cell_center, vel, index, m_xl, m_yl, m_xu, m_yu); }
  Grid upwindSurface(const Grid &cell_center, const Grid &vel, const int index, int xl, int yl, int xu, int yu);

  //Compute divergence of scalar quantity times velocity using Barton's method
  //Meant for advection terms
  Grid transportDivergence2D(const Grid &quantity, const std::vector<Grid> &vel) { return transportDivergence2D(quantity,vel,m_xl,m_yl,m_xu,m_yu); }
  Grid transportDivergence2D(const Grid &quantity, const std::vector<Grid> &vel, int xl, int yl, int xu, int yu);
  //Compute 1D derivative using Barton's method, for transport terms
  Grid transportDerivative1D(const Grid &quantity, const Grid &vel, const int index) { return transportDerivative1D(quantity,vel,index,m_xl,m_yl,m_xu,m_yu); }
  Grid transportDerivative1D(const Grid &quantity, const Grid &vel, const int index, int xl, int yl, int xu, int yu);

  //Compute divergence of vector quantity a_x,a_y
  Grid divergence2D(const Grid& a_x, const Grid& a_y){ return divergence2D(a_x, a_y, m_xl, m_yl, m_xu, m_yu); }
  Grid divergence2D(const Grid& a_x, const Grid& a_y, int xl, int yl, int xu, int yu);
  //Compute divergence of vector quantity a[0],a[1]
  Grid divergence2D(const std::vector<Grid>& a) { return divergence2D(a,m_xl,m_yl,m_xu,m_yu); }
  Grid divergence2D(const std::vector<Grid>& a, int xl, int yl, int xu, int yu);

  //Compute single-direction first derivative for non-transport term (central differencing)
  Grid derivative1D(const Grid &quantity, const int index){ return derivative1D(quantity,index, m_xl, m_yl, m_xu, m_yu); }
  Grid derivative1D(const Grid &quantity, const int index, int xl, int yl, int xu, int yu);

  //Compute single-direction first derivative for non-transport term (backward differencing)
  //When positive_forward is true, 'backward' is taken to mean 'in the direction of decreasing indices'
  //When false, 'backward' is taken to mean 'in the direction of increasing indices'
  Grid derivative1DBackward(const Grid &quantity, bool positive_forward, const int index){ return derivative1DBackward(quantity,positive_forward,index, m_xl, m_yl, m_xu, m_yu); }
  Grid derivative1DBackward(const Grid &quantity, bool positive_forward, const int index, int xl, int yl, int xu, int yu);

  //Compute single-direction second derivative
  Grid secondDerivative1D(const Grid &quantity, const int index) { return secondDerivative1D(quantity, index, m_xl, m_yl, m_xu, m_yu); }
  Grid secondDerivative1D(const Grid &quantity, const int index, int xl, int yl, int xu, int yu);

  Grid secondDerivative1DBackward(const Grid &quantity, bool positive_forward, const int index) { return secondDerivative1DBackward(quantity,positive_forward,index, m_xl, m_yl, m_xu, m_yu); }
  Grid secondDerivative1DBackward(const Grid &quantity, bool positive_forward, const int index, int xl, int yl, int xu, int yu);

  //Computes Laplacian (del squared) of "quantity"
  Grid laplacian(const Grid &quantity) { return laplacian(quantity,m_xl,m_yl,m_xu,m_yu); }
  Grid laplacian(const Grid &quantity, int xl, int yl, int xu, int yu);

  //Computes curl of vector in z-direction (result in xy-plane)
  std::vector<Grid> curlZ(const Grid& z) { return curlZ(z,m_xl,m_yl,m_xu,m_yu); }
  std::vector<Grid> curlZ(const Grid& z, int xl, int yl, int xu, int yu);

  //Computes curl of vector in xy-plane (result in z-direction)
  Grid curl2D(const Grid& x, const Grid& y) { return curl2D(x,y,m_xl,m_yl,m_xu,m_yu); }
  Grid curl2D(const Grid& x, const Grid& y, int xl, int yl, int xu, int yu);
  
  // boundary interpolation
  double boundaryInterpolate(const Grid &quantity, int i1, int j1, int i2, int j2);
  double boundaryExtrapolate(const Grid &quantity, int i1, int j1, int i2, int j2);

  // Savitzky-Golay differentiation
  Grid derivativeSGx(const Grid& grid) const {return m_sg.derivative1D(grid,0,m_grids[d_x].min());};
  Grid derivativeSGy(const Grid& grid) const {return m_sg.derivative1D(grid,1,m_grids[d_y].min());};

  template <typename T> static std::string num2str(T num,int precision=-1)
    {
        std::ostringstream ss;
        if(precision == -1) ss.precision(std::numeric_limits<double>::digits10 + 1);
        else ss.precision(precision);
        ss.precision(precision);
        ss << num;
        return ss.str();
    };
};

#endif