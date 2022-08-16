//evolution.cpp
//PlasmaDomain functionality relating to time-evolution of the plasma

#include "plasmadomain.hpp"
#include "constants.hpp"
#include <cmath>
// #include "solarutils.hpp"

void PlasmaDomain::run(double time_duration)
{
  // handle when time is supplied as input
  if(time_duration >= 0.0) m_duration = time_duration;
  else {
    assert(m_duration >= 0.0 && "Duration must be specified on command line, or in state file");
  }
  m_max_time = m_time + m_duration;
  
  // run loop for advanceTime
  while (m_time < m_max_time && (max_iterations < 0 || m_iter < max_iterations)){
    // advance time and note the time iteration before and after stepping
    int old_time_iter = (int)(m_time/m_time_output_interval);
    advanceTime();
    int new_time_iter = (int)(m_time/m_time_output_interval);
    // store data for output flags
    bool store_cond_1 = m_iter_output_interval > 0 && m_iter%m_iter_output_interval == 0;
    bool store_cond_2 = m_time_output_interval > 0.0 && new_time_iter > old_time_iter;
    if (store_cond_1 || store_cond_2) storeGrids();
    // write grid data to file
    if (m_output_counter%m_write_interval == 0){
      writeToOutFile();
      writeStateFile("end");
    }
  }  
  writeToOutFile();
  writeStateFile("end");
}

void PlasmaDomain::advanceTime(bool verbose)
{
  // determine timestep for this iteration
  double min_dt;
  double dt_raw = m_eqs->getDT().min(m_xl,m_yl,m_xu,m_yu);
  min_dt = epsilon*dt_raw;

  // iterate module functions
  m_module_handler.preIterateModules(min_dt);
  m_module_handler.iterateModules(min_dt);

  // determine viscosity coefficient
  double visc_coeff = epsilon_viscous*0.5*((m_grids[d_x].square() + m_grids[d_y].square())/dt_raw).min();

  // run time integration according to specified method
  if (m_time_integrator == TimeIntegrator::RK2) integrateRK2(min_dt, visc_coeff);
  else if (m_time_integrator == TimeIntegrator::RK4) integrateRK4(min_dt, visc_coeff);
  else if (m_time_integrator == TimeIntegrator::Euler) integrateEuler(min_dt, visc_coeff);

  // post iterate all modules
  m_module_handler.postIterateModules(min_dt);

  // print update to terminal
  if(m_std_out_interval > 0 && m_iter%m_std_out_interval == 0) printUpdate(min_dt);

  // step time
  m_time += min_dt;
  m_iter++;
}

void PlasmaDomain::integrateEuler(double time_step, double visc_coeff)
{
  std::vector<Grid> time_derivatives = m_eqs->computeTimeDerivatives(visc_coeff);
  m_eqs->applyTimeDerivatives(time_derivatives, time_step);
}

void PlasmaDomain::integrateRK2(double time_step, double visc_coeff)
{
  // First create copy of system advanced by half the timestep (Euler)...
  std::vector<Grid> time_derivatives_naive = m_eqs->computeTimeDerivatives(visc_coeff);
  std::vector<Grid> grids_halfstep = m_eqs->allGrids();
  m_eqs->applyTimeDerivatives(grids_halfstep, time_derivatives_naive, 0.5*time_step);

  // ...then use that half-advanced system to compute time derivatives to
  // apply to actual system for full time step (second-order R-K, a.k.a. midpoint method)
  std::vector<Grid> time_derivatives = m_eqs->computeTimeDerivatives(grids_halfstep,visc_coeff);
  m_eqs->applyTimeDerivatives(time_derivatives, time_step);
}

void PlasmaDomain::integrateRK4(double time_step, double visc_coeff)
{
  std::vector<Grid> k1 = m_eqs->computeTimeDerivatives(visc_coeff);

  std::vector<Grid> grids_copy = m_eqs->allGrids();
  m_eqs->applyTimeDerivatives(grids_copy, k1, 0.5*time_step);
  std::vector<Grid> k2 = m_eqs->computeTimeDerivatives(grids_copy,visc_coeff);

  grids_copy = m_eqs->allGrids();
  m_eqs->applyTimeDerivatives(grids_copy, k2, 0.5*time_step);
  std::vector<Grid> k3 = m_eqs->computeTimeDerivatives(grids_copy,visc_coeff);

  grids_copy = m_eqs->allGrids();
  m_eqs->applyTimeDerivatives(grids_copy, k3, time_step);
  std::vector<Grid> k4 = m_eqs->computeTimeDerivatives(grids_copy,visc_coeff);

  std::vector<Grid> k_final = k1;
  for(int i=0; i<k_final.size(); i++)
    k_final[i] = (k1[i] + k4[i])/6.0 + (k2[i] + k3[i])/3.0;
  
  m_eqs->applyTimeDerivatives(k_final, time_step);
}

void PlasmaDomain::updateGhostZones()
{
  if(x_bound_1 != BoundaryCondition::Periodic && x_bound_1 != BoundaryCondition::OpenMoC){
    if(x_bound_1 == BoundaryCondition::Open) for(int j=m_yl; j<=m_yu; j++) openBoundaryExtrapolate(0, 1, 2, 3, j, j, j, j);
    else if(x_bound_1 == BoundaryCondition::Reflect) for(int j=m_yl; j<=m_yu; j++) reflectBoundaryExtrapolate(0, 1, 2, 3, j, j, j, j);
    else if(x_bound_1 == BoundaryCondition::Fixed) for(int j=m_yl; j<=m_yu; j++) fixedBoundaryExtrapolate(0, 1, 2, 3, j, j, j, j);
  }
  if(x_bound_2 != BoundaryCondition::Periodic && x_bound_2 != BoundaryCondition::OpenMoC){
    if(x_bound_2 == BoundaryCondition::Open) for(int j=m_yl; j<=m_yu; j++) openBoundaryExtrapolate(m_xdim-1, m_xdim-2, m_xdim-3, m_xdim-4, j, j, j, j);
    else if(x_bound_2 == BoundaryCondition::Reflect) for(int j=m_yl; j<=m_yu; j++) reflectBoundaryExtrapolate(m_xdim-1, m_xdim-2, m_xdim-3, m_xdim-4, j, j, j, j);
    else if(x_bound_2 == BoundaryCondition::Fixed) for(int j=m_yl; j<=m_yu; j++) fixedBoundaryExtrapolate(m_xdim-1, m_xdim-2, m_xdim-3, m_xdim-4, j, j, j, j);
  }
  if(y_bound_1 != BoundaryCondition::Periodic && y_bound_1 != BoundaryCondition::OpenMoC){
    if(y_bound_1 == BoundaryCondition::Open) for(int i=m_xl; i<=m_xu; i++) openBoundaryExtrapolate(i, i, i, i, 0, 1, 2, 3);
    else if(y_bound_1 == BoundaryCondition::Reflect) for(int i=m_xl; i<=m_xu; i++) reflectBoundaryExtrapolate(i, i, i, i, 0, 1, 2, 3);
    else if(y_bound_1 == BoundaryCondition::Fixed) for(int i=m_xl; i<=m_xu; i++) fixedBoundaryExtrapolate(i, i, i, i, 0, 1, 2, 3);
  }
  if(y_bound_2 != BoundaryCondition::Periodic && y_bound_2 != BoundaryCondition::OpenMoC){
    if(y_bound_2 == BoundaryCondition::Open) for(int i=m_xl; i<=m_xu; i++) openBoundaryExtrapolate(i, i, i, i, m_ydim-1, m_ydim-2, m_ydim-3, m_ydim-4);
    else if(y_bound_2 == BoundaryCondition::Reflect) for(int i=m_xl; i<=m_xu; i++) reflectBoundaryExtrapolate(i, i, i, i, m_ydim-1, m_ydim-2, m_ydim-3, m_ydim-4);
    else if(y_bound_2 == BoundaryCondition::Fixed) for(int i=m_xl; i<=m_xu; i++) fixedBoundaryExtrapolate(i, i, i, i, m_ydim-1, m_ydim-2, m_ydim-3, m_ydim-4);
  }
}

//Applies the open boundary condition to the cells indexed by the given indices
//assuming that (i1,j1) is at the edge of the domain and (i4,j4) is on the interior
//Assumes two ghost zones (i.e. N_GHOST == 2), will abort if not
//Meant to be called repeatedly during advanceTime; energy is handled, not temperature
void PlasmaDomain::openBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4)
{
  assert(N_GHOST == 2 && "This function assumes two ghost zones");
  assert((i1 == i2 || j1 == j2) && "Points to extrapolate must be x-aligned or y-aligned");
  bool x_boundary = (j1 == j2);

  // Grid &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y], &m_rho = grids[rho], &m_thermal_energy = grids[thermal_energy];

  double delta_last = x_boundary ? m_grids[d_x](i3,j3) : m_grids[d_y](i3,j3);
  double dist23 = x_boundary ? 0.5*(m_grids[d_x](i2,j2) + m_grids[d_x](i3,j3)) : 0.5*(m_grids[d_y](i2,j2) + m_grids[d_y](i3,j3));
  double dist12 = x_boundary ? 0.5*(m_grids[d_x](i1,j1) + m_grids[d_x](i2,j2)) : 0.5*(m_grids[d_y](i1,j1) + m_grids[d_y](i2,j2));
  double scale_2 = std::pow(open_boundary_decay_base,dist23/delta_last);
  double scale_1 = std::pow(open_boundary_decay_base,dist12/delta_last);
  for(int v : m_eqs->densities()){
    Grid& rho = m_eqs->grid(v);
    rho(i1,j1) = scale_1*rho(i3,j3);
    rho(i2,j2) = scale_2*rho(i3,j3);
  }
  for(int v : m_eqs->thermal_energies()){
    Grid& energy = m_eqs->grid(v);
    energy(i1,j1) = scale_1*energy(i3,j3);
    energy(i2,j2) = scale_2*energy(i3,j3);
  }

  // get maximum c_s
  double c_s{0.};
  for(int species = 0; species < m_eqs->num_species(); species++){
    Grid m_press = (m_adiabatic_index - 1.0)*m_eqs->grid(m_eqs->thermal_energies()[species]);
    Grid& m_rho = m_eqs->grid(m_eqs->densities()[species]);
    double c_s_new = std::sqrt(m_adiabatic_index*m_press(i3,j3)/m_rho(i3,j3));
    if (c_s_new > c_s) c_s = c_s_new;
  }

  for(int species = 0; species < m_eqs->num_species(); species++){
    Grid& m_mom_x = m_eqs->grid(m_eqs->momenta()[species][0]);
    Grid& m_mom_y = m_eqs->grid(m_eqs->momenta()[species][1]);
    Grid& m_rho = m_eqs->grid(m_eqs->densities()[species]);

    double vel_x = m_mom_x(i3,j3)/m_rho(i3,j3);
    double vel_y = m_mom_y(i3,j3)/m_rho(i3,j3);

    Grid m_press = (m_adiabatic_index - 1.0)*m_eqs->grid(m_eqs->thermal_energies()[species]);

    double boost_vel = open_boundary_strength*c_s;
    if(i2 > i1 || j2 > j1) boost_vel *= -1.0;
    
    //Determine ghost cell velocity s.t. interpolated velocity at interior cell boundary is some multiple of sound speed, outward
    double dist = dist23;
    if(x_boundary){
      double boundary_vel;
      if(i1 > i2) boundary_vel = std::max(0.0, vel_x + boost_vel);
      else { assert(i2 > i1); boundary_vel = std::min(0.0, vel_x + boost_vel); }
      double ghost_vel = (dist*(boundary_vel) - 0.5*m_grids[d_x](i2,j2)*vel_x)/(0.5*m_grids[d_x](i3,j3));
      m_mom_x(i1,j1) = m_rho(i1,j1)*ghost_vel;
      m_mom_x(i2,j2) = m_rho(i2,j2)*ghost_vel;
      m_mom_y(i1,j1) = m_rho(i1,j1)*vel_y;
      m_mom_y(i2,j2) = m_rho(i2,j2)*vel_y;
    } else {
      assert(i1 == i2);
      double boundary_vel;
      if(j1 > j2) boundary_vel = std::max(0.0, vel_y + boost_vel);
      else { assert(j2 > j1); boundary_vel = std::min(0.0, vel_y + boost_vel); }
      double ghost_vel = (dist*(boundary_vel) - 0.5*m_grids[d_y](i2,j2)*vel_y)/(0.5*m_grids[d_y](i3,j3)); //add vel_y to c_s? sound speed relative to current bulk velocity?
      m_mom_x(i1,j1) = m_rho(i1,j1)*vel_x;
      m_mom_x(i2,j2) = m_rho(i2,j2)*vel_x;
      m_mom_y(i1,j1) = m_rho(i1,j1)*ghost_vel;
      m_mom_y(i2,j2) = m_rho(i2,j2)*ghost_vel;
    }
  }
}

//Applies the reflect boundary condition to the cells indexed by the given indices
//assuming that (i1,j1) is at the edge of the domain and (i4,j4) is on the interior
//Matches the thermal energy and mass density of the closest interior (i.e. non-ghost) cell
//throughout the ghost zone and halts advection at the boundary
//Assumes two ghost zones (i.e. N_GHOST == 2), will abort if not
void PlasmaDomain::reflectBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4)
{
  assert(N_GHOST == 2 && "This function assumes two ghost zones");
  // Grid &m_rho = grids[rho], &m_thermal_energy = grids[thermal_energy], &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y];
  // m_rho(i1,j1) = m_rho(i3,j3); m_rho(i2,j2) = m_rho(i3,j3); //Match density of nearest interior cell
  // m_thermal_energy(i1,j1) = m_thermal_energy(i3,j3); m_thermal_energy(i2,j2) = m_thermal_energy(i3,j3); //Match temperature of nearest interior cell
  for(int v : m_eqs->thermal_energies()){
    m_eqs->grid(v)(i1,j1) = m_eqs->grid(v)(i3,j3);
    m_eqs->grid(v)(i2,j2) = m_eqs->grid(v)(i3,j3); //Match thermal energy of nearest interior cell
  }
  for(int v : m_eqs->densities()){
    m_eqs->grid(v)(i1,j1) = m_eqs->grid(v)(i3,j3);
    m_eqs->grid(v)(i2,j2) = m_eqs->grid(v)(i3,j3); //Match density of nearest interior cell
  }
  // m_mom_x(i1,j1) = 0.0; m_mom_x(i2,j2) = 0.0; m_mom_x(i3,j3) = 0.0;
  // m_mom_y(i1,j1) = 0.0; m_mom_y(i2,j2) = 0.0; m_mom_y(i3,j3) = 0.0;
  for(std::vector<int> mom : m_eqs->momenta()){
    for(int v : mom){
      m_eqs->grid(v)(i1,j1) = 0.0;
      m_eqs->grid(v)(i2,j2) = 0.0;
      m_eqs->grid(v)(i3,j3) = 0.0;
    }
  }
}

//Applies the fixed boundary condition to the cells indexed by the given indices
//assuming that (i1,j1) is at the edge of the domain and (i4,j4) is on the interior
//Halts all advection at the boundary; leaves all other quantities in the ghost zones fixed
//Assumes two ghost zones (i.e. N_GHOST == 2), will abort if not
void PlasmaDomain::fixedBoundaryExtrapolate(int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4)
{
  assert(N_GHOST == 2 && "This function assumes two ghost zones");
  // Grid &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y];
  // m_mom_x(i1,j1) = 0.0; m_mom_x(i2,j2) = 0.0; m_mom_x(i3,j3) = 0.0;
  // m_mom_y(i1,j1) = 0.0; m_mom_y(i2,j2) = 0.0; m_mom_y(i3,j3) = 0.0;
  for(std::vector<int> mom : m_eqs->momenta()){
    for(int v : mom){
      m_eqs->grid(v)(i1,j1) = 0.0;
      m_eqs->grid(v)(i2,j2) = 0.0;
      m_eqs->grid(v)(i3,j3) = 0.0;
    }
  }
}