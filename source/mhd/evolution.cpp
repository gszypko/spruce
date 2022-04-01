//evolution.cpp
//PlasmaDomain functionality relating to time-evolution of the plasma

#include "plasmadomain.hpp"
#include "solarutils.hpp"

void PlasmaDomain::run(double time_duration)
{
  max_time = m_time + time_duration;
  if(!continue_mode){
    outputPreamble();
    outputCurrentState();
  }
  while( m_time < max_time && (max_iterations < 0 || m_iter < max_iterations)){
    #if BENCHMARKING_ON
    InstrumentationTimer timer((std::string("iteration ") + std::to_string(m_iter)).c_str());
    #endif
    double old_time = m_time;
    advanceTime();
    if(sg_filter_interval > 0 && m_iter != 0 && m_iter%sg_filter_interval == 0) filterSavitzkyGolay();
    if((iter_output_interval > 0 && m_iter%iter_output_interval == 0) || (time_output_interval > 0.0
        && (int)(m_time/time_output_interval) > (int)(old_time/time_output_interval))) outputCurrentState();
    if(safe_state_mode) writeStateFile();
  }
  if(safe_state_mode) cleanUpStateFiles();
  else writeStateFile();
}

void PlasmaDomain::advanceTime(bool verbose)
{
  double min_dt;

  double dt_raw = m_grids[dt].min(m_xl,m_yl,m_xu,m_yu);
  min_dt = epsilon*dt_raw;

  m_module_handler.iterateModules(min_dt);

  double visc_coeff = epsilon_viscous*0.5*((m_grids[d_x].square() + m_grids[d_y].square())/dt_raw).min();

  if(time_integrator == TimeIntegrator::RK2) integrateRK2(m_grids, min_dt, visc_coeff);
  else if(time_integrator == TimeIntegrator::RK4) integrateRK4(m_grids, min_dt, visc_coeff);
  else if(time_integrator == TimeIntegrator::Euler) integrateEuler(m_grids, min_dt, visc_coeff);

  propagateChanges();

  if(std_out_interval > 0 && m_iter%std_out_interval == 0) printUpdate(min_dt);

  m_time += min_dt;
  m_iter++;
}

void PlasmaDomain::integrateEuler(std::vector<Grid> &grids, double time_step, double visc_coeff)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  std::vector<Grid> time_derivatives = computeTimeDerivatives(grids,visc_coeff);
  applyTimeDerivatives(grids, time_derivatives, time_step);
}

void PlasmaDomain::integrateRK2(std::vector<Grid> &grids, double time_step, double visc_coeff)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  // First create copy of system advanced by half the timestep (Euler)...
  std::vector<Grid> time_derivatives_naive = computeTimeDerivatives(grids,visc_coeff);
  std::vector<Grid> grids_halfstep = grids;
  applyTimeDerivatives(grids_halfstep, time_derivatives_naive, 0.5*time_step);

  // ...then use that half-advanced system to compute time derivatives to
  // apply to actual system for full time step (second-order R-K, a.k.a. midpoint method)
  std::vector<Grid> time_derivatives = computeTimeDerivatives(grids_halfstep,visc_coeff);
  applyTimeDerivatives(grids, time_derivatives, time_step);
}

void PlasmaDomain::integrateRK4(std::vector<Grid> &grids, double time_step, double visc_coeff)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  std::vector<Grid> k1 = computeTimeDerivatives(grids,visc_coeff);

  std::vector<Grid> grids_copy = grids;
  applyTimeDerivatives(grids_copy, k1, 0.5*time_step);
  std::vector<Grid> k2 = computeTimeDerivatives(grids_copy,visc_coeff);

  grids_copy = grids;
  applyTimeDerivatives(grids_copy, k2, 0.5*time_step);
  std::vector<Grid> k3 = computeTimeDerivatives(grids_copy,visc_coeff);

  grids_copy = grids;
  applyTimeDerivatives(grids_copy, k3, time_step);
  std::vector<Grid> k4 = computeTimeDerivatives(grids_copy,visc_coeff);

  std::vector<Grid> k_final = k1;
  for(int i=0; i<k_final.size(); i++)
    k_final[i] = (k1[i] + k4[i])/6.0 + (k2[i] + k3[i])/3.0;
  
  applyTimeDerivatives(grids, k_final, time_step);
}

// Adds time derivatives (in order of {d_rho_dt, d_mom_x_dt, d_mom_y_dt, d_bi_x_dt, d_bi_y_dt, d_thermal_energy_dt})
// to the Grids of the system given (in place), multiplied by the given step
// and applies all post-step cleanup (floor underdensity, recompute derived quantities, update ghost zones)
void PlasmaDomain::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  grids[rho] += step*m_ghost_zone_mask*time_derivatives[0];
  grids[mom_x] += step*m_ghost_zone_mask*time_derivatives[1];
  grids[mom_y] += step*m_ghost_zone_mask*time_derivatives[2];
  grids[bi_x] += step*m_ghost_zone_mask*time_derivatives[3];
  grids[bi_y] += step*m_ghost_zone_mask*time_derivatives[4];
  grids[thermal_energy] += step*m_ghost_zone_mask*time_derivatives[5];
  propagateChanges(grids);
}

// Returns time derivatives in the following order:
// {d_rho_dt, d_mom_x_dt, d_mom_y_dt, d_bi_x_dt, d_bi_y_dt, d_thermal_energy_dt}
std::vector<Grid> PlasmaDomain::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  //Advance time by min_dt
  const Grid &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y],
        &m_v_x = grids[v_x], &m_v_y = grids[v_y],
        &m_be_x = grids[be_x], &m_be_y = grids[be_y],
        &m_bi_x = grids[bi_x], &m_bi_y = grids[bi_y],
        &m_rho = grids[rho], &m_thermal_energy = grids[thermal_energy], &m_press = grids[press];

  std::vector<Grid> m_v = {m_v_x,m_v_y};

  Grid viscous_force_x = visc_coeff*laplacian(m_mom_x);
  Grid viscous_force_y = visc_coeff*laplacian(m_mom_y);

  Grid d_rho_dt = -transportDivergence2D(m_rho,m_v);

  Grid curl_db = curl2D(m_bi_x,m_bi_y)/(4.0*PI);
  std::vector<Grid> mag_mom_terms_firstorder = Grid::CrossProductZ2D(curl_db,{m_be_x,m_be_y});
  std::vector<Grid> mag_mom_terms_secorder = Grid::CrossProductZ2D(curl_db,{m_bi_x,m_bi_y}); //second order terms
  Grid d_mom_x_dt = -transportDivergence2D(m_mom_x, m_v)
                  - derivative1D(m_press, 0)
                  + m_ghost_zone_mask * (m_rho*grids[grav_x] + viscous_force_x
                  + mag_mom_terms_firstorder[0] + mag_mom_terms_secorder[0]);
  Grid d_mom_y_dt = -transportDivergence2D(m_mom_y, m_v)
                  - derivative1D(m_press, 1)
                  + m_ghost_zone_mask * (m_rho*grids[grav_y] + viscous_force_y
                  + mag_mom_terms_firstorder[1] + mag_mom_terms_secorder[1]);

  std::vector<Grid> induction_rhs_b0 = curlZ(Grid::CrossProduct2D({m_v_x,m_v_y},{m_be_x,m_be_y}));
  std::vector<Grid> induction_rhs_db = curlZ(Grid::CrossProduct2D({m_v_x,m_v_y},{m_bi_x,m_bi_y}));
  Grid d_bi_x_dt = induction_rhs_b0[0] + induction_rhs_db[0];
  Grid d_bi_y_dt = induction_rhs_b0[1] + induction_rhs_db[1];

  // // Grid mag_energy_term = computeMagneticEnergyTerm();
  // // Grid d_thermal_energy_dt = - transportDivergence2D(m_thermal_energy, m_v)
  // //                          - transportDivergence2D(m_grids[kinetic_energy], m_v)
  // //                          - transportDivergence2D(m_grids[mag_energy], m_v)
  // //                          - divergence2D(m_press*m_v_x, m_press*m_v_y)
  // //                          + (m_rho*m_grids[grav_x] + viscous_force_x)*m_v_x
  // //                          + (m_rho*m_grids[grav_y] + viscous_force_y)*m_v_y
  // //                          + 0.5*(m_v_x.square() + m_v_y.square())*d_rho_dt
  // //                          - m_v_x*d_mom_x_dt - m_v_y*d_mom_y_dt
  // //                          - Grid::DotProduct2D({m_b_x+m_db_x,m_b_y+m_db_y},{d_db_x_dt,d_db_y_dt})/(4.0*PI)
  // //                          + mag_energy_term;
  // Grid d_thermal_energy_dt = - transportDivergence2D(m_thermal_energy,m_v)
  //                            - m_press*divergence2D(m_v);

  // // //NOTE: conduction, radiation, heating should also be applied to this equation
  // // double nu_i_e; //electron-ion collision frequency
  // // Grid m_thermal_energy_e;
  // // Grid m_press_e;
  // // Grid temp_e;
  // // Grid m_thermal_energy_i;
  // // Grid m_press_i;
  // // Grid temp_i;

  // // Grid collision_term = K_B*m_grids[DerivedVars::n]/(gamma - 1.0)*nu_i_e*(temp_i - temp_e);

  // // //Can still add: electric potential term
  // // Grid d_thermal_energy_e_dt = - transportDerivative1D(m_thermal_energy_e, m_v_x, 0)
  // //                          - transportDerivative1D(m_thermal_energy_e, m_v_y, 1)
  // //                          - derivative1D(m_press_e*m_v_x, 0) - derivative1D(m_press_e*m_v_y, 1)
  // //                          + (viscous_force_x)*m_v_x
  // //                          + (viscous_force_y)*m_v_y
  // //                          + collision_term;
  
  // // //NOTE: heating may or may not be applied to this equation
  // // //Can still add: spitzer viscosity, electric potential term
  // // Grid d_thermal_energy_i_dt = - transportDerivative1D(m_thermal_energy_i + m_grids[kinetic_energy], m_v_x, 0)
  // //                          - transportDerivative1D(m_thermal_energy_i + m_grids[kinetic_energy], m_v_y, 1)
  // //                          - derivative1D(m_press_i*m_v_x, 0) - derivative1D(m_press_i*m_v_y, 1)
  // //                          - derivative1D(m_grids[mag_pxx]*m_v_x + m_grids[mag_pxy]*m_v_y, 0)
  // //                          - derivative1D(m_grids[mag_pxy]*m_v_x + m_grids[mag_pyy]*m_v_y, 1)
  // //                          + (m_rho*m_grids[grav_x] + viscous_force_x)*m_v_x
  // //                          + (m_rho*m_grids[grav_y] + viscous_force_y)*m_v_y
  // //                          + 0.5*(m_v_x.square() + m_v_y.square())*d_rho_dt
  // //                          - m_v_x*d_mom_x_dt - m_v_y*d_mom_y_dt
  // //                          - collision_term;

  Grid d_thermal_energy_dt = - transportDivergence2D(m_thermal_energy,m_v)
                             - m_press*divergence2D(m_v);
  // if(ambient_heating) d_thermal_energy_dt += m_ghost_zone_mask*heating_rate;

  return {d_rho_dt, d_mom_x_dt, d_mom_y_dt, d_bi_x_dt, d_bi_y_dt, d_thermal_energy_dt};
}

Grid PlasmaDomain::computeMagneticEnergyTerm()
{
  Grid &m_b_x = m_grids[be_x], &m_b_y = m_grids[be_y],
       &m_db_x = m_grids[bi_x], &m_db_y = m_grids[bi_y],
       &m_v_x = m_grids[v_x], &m_v_y = m_grids[v_y];
  std::vector<Grid> m_b = {m_b_x,m_b_y}, m_db = {m_db_x,m_db_y}, m_v = {m_v_x,m_v_y};
  Grid b_sq = Grid::DotProduct2D(m_b,m_b), db_sq = Grid::DotProduct2D(m_db,m_db),
       b_dot_v = Grid::DotProduct2D(m_b,m_v), db_dot_v = Grid::DotProduct2D(m_db,m_v),
       b_dot_db = Grid::DotProduct2D(m_b,m_db);

  std::vector<Grid> bg_field_flux = {0.5*m_v_x*b_sq - m_b_x*b_dot_v, 0.5*m_v_y*b_sq - m_b_y*b_dot_v};
  std::vector<Grid> db_field_flux = {0.5*m_v_x*db_sq - m_db_x*db_dot_v, 0.5*m_v_y*db_sq - m_db_y*db_dot_v}; //second order
  std::vector<Grid> crossterm_flux = {m_v_x*b_dot_db - m_b_x*db_dot_v - m_db_x*b_dot_v, m_v_y*b_dot_db - m_b_y*db_dot_v - m_db_y*b_dot_v};

  return -1.0/(4.0*PI)*(divergence2D(bg_field_flux) + divergence2D(db_field_flux) + divergence2D(crossterm_flux));

}

//Computes the magnetic magnitude, direction, and pressure terms
//based on the current values of be_x, be_y
//Also computes grid spacing d_x and d_y from pos_x and pos_y
//For terms that are computed once, at initialization, and unchanged thereafter
void PlasmaDomain::computeConstantTerms()
{
  m_grids[div_be] = divergence2D(m_grids[be_x],m_grids[be_y]);
}

//Propagates all changes to mass, momentum, energy, and magnetic field
//to all other quantities s.t. everything is consistent
//Includes updating ghost zones
void PlasmaDomain::propagateChanges()
{
  propagateChanges(m_grids);
}

void PlasmaDomain::propagateChanges(std::vector<Grid> &grids)
{
  catchUnderdensity(grids);
  recomputeTemperature(grids);
  recomputeDerivedVariables(grids);
  updateGhostZones(grids);
}

//Recompute pressure, energy, velocity, radiative loss rate, dt, dt_thermal, dt_rad
//from the current base variables. Does not compute variables that are shut off
//by the member physics settings.
void PlasmaDomain::recomputeDerivedVariables()
{
  recomputeDerivedVariables(m_grids);
}

void PlasmaDomain::recomputeDerivedVariables(std::vector<Grid> &grids)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  grids[press] = 2.0*K_B*grids[rho]*grids[temp]/m_ion_mass;
  grids[thermal_energy] = grids[press]/(m_adiabatic_index - 1.0);
  grids[kinetic_energy] = 0.5*(grids[mom_x].square() + grids[mom_y].square())/grids[rho];
  grids[v_x] = grids[mom_x]/grids[rho];
  grids[v_y] = grids[mom_y]/grids[rho];
  grids[n] = grids[rho]/m_ion_mass;
  grids[div_bi] = divergence2D(grids[bi_x],grids[bi_y]);
  Grid m_b_x = grids[be_x] + grids[bi_x], m_b_y = grids[be_y] + grids[bi_y];
  grids[b_magnitude] = (m_b_x.square() + m_b_y.square()).sqrt();
  //Need to ensure that b_hat is zero when b is zero
  Grid &b_h_x = grids[b_hat_x], &b_h_y = grids[b_hat_y], &b_mag = grids[b_magnitude];
  b_h_x = m_b_x/b_mag;
  b_h_y = m_b_y/b_mag;
  for(int i=0; i<m_xdim; i++){
    for(int j=0; j<m_ydim; j++){
      if(b_mag(i,j) == 0.0){
        b_h_x(i,j) = 0.0;
        b_h_y(i,j) = 0.0;
      }
    }
  }
  recomputeDT();
}

void PlasmaDomain::recomputeTemperature()
{
  recomputeTemperature(m_grids);
}

void PlasmaDomain::recomputeTemperature(std::vector<Grid> &grids)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  Grid &m_press = grids[press], &m_thermal_energy = grids[thermal_energy];
  m_thermal_energy = m_thermal_energy.max(thermal_energy_min);
  m_press = (m_adiabatic_index - 1.0)*m_thermal_energy;
  grids[temp] = (m_ion_mass*m_press/(2.0*K_B*grids[rho])).max(temp_min);
}

void PlasmaDomain::catchUnderdensity()
{
  catchUnderdensity(m_grids);
}

void PlasmaDomain::catchUnderdensity(std::vector<Grid> &grids)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  for(int i=0; i<m_xdim; i++){
    for(int j=0; j<m_ydim; j++){
      if(grids[rho](i,j) < rho_min){
        //Only notify if not within ghost zones
        if(i >= m_xl && i <= m_xu && j >= m_yl && j <= m_yu){
          // std::cout << "Density too low at (" << i << "," << j << ")" << "\n";
          grids[rho](i,j) = rho_min;
        }
      }
    }
  }
}

//Try taking diagonal of grid cell
//and then magnitude of velocity
//with safety coefficient
void PlasmaDomain::recomputeDT()
{
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid c_s = (m_adiabatic_index*m_grids[press]/m_grids[rho]).sqrt();
  Grid c_s_sq = c_s.square();

  Grid v_alfven = m_grids[b_magnitude]/(4.0*PI*m_grids[rho]).sqrt();
  m_grids[v_a] = v_alfven;
  Grid v_alfven_sq = v_alfven.square();
  Grid one = Grid::Ones(m_xdim,m_ydim);

  //Magnetoacoustic modes
  Grid delta = (one - 4.0*c_s_sq*v_alfven_sq/(c_s_sq+v_alfven_sq).square()).sqrt();
  Grid v_fast = (0.5*(c_s_sq + v_alfven_sq)*(one + delta)).sqrt();
  Grid v_slow = (0.5*(c_s_sq + v_alfven_sq)*(one - delta)).sqrt();

  //Bulk velocity transit time
  Grid diagonals = (m_grids[d_x].square() + m_grids[d_y].square()).sqrt();
  Grid vel_mag = (m_grids[v_x].square() + m_grids[v_y].square()).sqrt();

  m_grids[dt] = diagonals/(c_s + v_alfven + v_fast + v_slow + vel_mag);
}


void PlasmaDomain::updateGhostZones()
{
  updateGhostZones(m_grids);
}

void PlasmaDomain::updateGhostZones(std::vector<Grid> &grids)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  if(x_bound_1 != BoundaryCondition::Periodic){
    for(int j=m_yl; j<=m_yu; j++){
      if(x_bound_1 == BoundaryCondition::Open) openBoundaryExtrapolate(grids, 0, 1, 2, 3, j, j, j, j);
      else if(x_bound_1 == BoundaryCondition::Reflect) reflectBoundaryExtrapolate(grids, 0, 1, 2, 3, j, j, j, j);
      else if(x_bound_1 == BoundaryCondition::Fixed) fixedBoundaryExtrapolate(grids, 0, 1, 2, 3, j, j, j, j);
    }
  }
  if(x_bound_2 != BoundaryCondition::Periodic){
    for(int j=m_yl; j<=m_yu; j++){
      if(x_bound_2 == BoundaryCondition::Open) openBoundaryExtrapolate(grids, m_xdim-1, m_xdim-2, m_xdim-3, m_xdim-4, j, j, j, j);
      else if(x_bound_2 == BoundaryCondition::Reflect) reflectBoundaryExtrapolate(grids, m_xdim-1, m_xdim-2, m_xdim-3, m_xdim-4, j, j, j, j);
      else if(x_bound_2 == BoundaryCondition::Fixed) fixedBoundaryExtrapolate(grids, m_xdim-1, m_xdim-2, m_xdim-3, m_xdim-4, j, j, j, j);
    }
  }
  if(y_bound_1 != BoundaryCondition::Periodic){
    for(int i=m_xl; i<=m_xu; i++){
      if(y_bound_1 == BoundaryCondition::Open) openBoundaryExtrapolate(grids, i, i, i, i, 0, 1, 2, 3);
      else if(y_bound_1 == BoundaryCondition::Reflect) reflectBoundaryExtrapolate(grids, i, i, i, i, 0, 1, 2, 3);
      else if(y_bound_1 == BoundaryCondition::Fixed) fixedBoundaryExtrapolate(grids, i, i, i, i, 0, 1, 2, 3);
    }
  }
  if(y_bound_2 != BoundaryCondition::Periodic){
    for(int i=m_xl; i<=m_xu; i++){
      if(y_bound_2 == BoundaryCondition::Open) openBoundaryExtrapolate(grids, i, i, i, i, m_ydim-1, m_ydim-2, m_ydim-3, m_ydim-4);
      else if(y_bound_2 == BoundaryCondition::Reflect) reflectBoundaryExtrapolate(grids, i, i, i, i, m_ydim-1, m_ydim-2, m_ydim-3, m_ydim-4);
      else if(y_bound_2 == BoundaryCondition::Fixed) fixedBoundaryExtrapolate(grids, i, i, i, i, m_ydim-1, m_ydim-2, m_ydim-3, m_ydim-4);
    }
  }
}

//Applies the open boundary condition to the cells indexed by the given indices
//assuming that (i1,j1) is at the edge of the domain and (i4,j4) is on the interior
//Assumes two ghost zones (i.e. N_GHOST == 2), will abort if not
//Meant to be called repeatedly during advanceTime; energy is handled, not temperature
void PlasmaDomain::openBoundaryExtrapolate(std::vector<Grid> &grids, int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  assert(N_GHOST == 2 && "This function assumes two ghost zones");
  assert((i1 == i2 || j1 == j2) && "Points to extrapolate must be x-aligned or y-aligned");
  Grid &m_pos_x = grids[pos_x], &m_pos_y = grids[pos_y];

  Grid &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y], &m_rho = grids[rho], &m_thermal_energy = grids[thermal_energy];

  // m_rho(i1,j1) = 0.25*m_rho(i3,j3); m_rho(i2,j2) = 0.5*m_rho(i3,j3);
  // m_thermal_energy(i1,j1) = 0.25*m_thermal_energy(i3,j3); m_thermal_energy(i2,j2) = 0.5*m_thermal_energy(i3,j3);
  double delta_last = std::abs(m_pos_x(i3,j3) - m_pos_x(i4,j4) + m_pos_y(i3,j3) - m_pos_y(i4,j4));
  double scale_2 = std::pow(open_boundary_decay_base,std::abs(m_pos_x(i2,j2) - m_pos_x(i3,j3) + m_pos_y(i2,j2) - m_pos_y(i3,j3))/delta_last);
  double scale_1 = std::pow(open_boundary_decay_base,std::abs(m_pos_x(i1,j1) - m_pos_x(i3,j3) + m_pos_y(i1,j1) - m_pos_y(i3,j3))/delta_last);
  m_rho(i1,j1) = scale_1*m_rho(i3,j3); m_rho(i2,j2) = scale_2*m_rho(i3,j3);
  m_thermal_energy(i1,j1) = scale_1*m_thermal_energy(i3,j3); m_thermal_energy(i2,j2) = scale_2*m_thermal_energy(i3,j3);

  double vel_x = m_mom_x(i3,j3)/m_rho(i3,j3);
  double vel_y = m_mom_y(i3,j3)/m_rho(i3,j3);

  // double boundary_strength = 1.0;
  double c_s = std::sqrt(m_adiabatic_index*grids[press](i3,j3)/m_rho(i3,j3));
  double boost_vel = open_boundary_strength*c_s;
  if(i2 > i1 || j2 > j1) boost_vel *= -1.0;
  
  //Determine ghost cell velocity s.t. interpolated velocity at interior cell boundary is some multiple of sound speed, outward
  double dist = std::abs(m_pos_x(i3,j3) - m_pos_x(i2,j2) + m_pos_y(i3,j3) - m_pos_y(i2,j2));
  if(i1 == i2){
    double boundary_vel;
    if(j1 > j2) boundary_vel = std::max(0.0, vel_y + boost_vel);
    else { assert(j2 > j1); boundary_vel = std::min(0.0, vel_y + boost_vel); }
    double ghost_vel = (dist*(boundary_vel) - 0.5*grids[d_y](i2,j2)*vel_y)/(0.5*grids[d_y](i3,j3)); //add vel_y to c_s? sound speed relative to current bulk velocity?
    m_mom_x(i1,j1) = m_rho(i1,j1)*vel_x;
    m_mom_x(i2,j2) = m_rho(i2,j2)*vel_x;
    m_mom_y(i1,j1) = m_rho(i1,j1)*ghost_vel;
    m_mom_y(i2,j2) = m_rho(i2,j2)*ghost_vel;
  } else {
    assert(j1 == j2);
    double boundary_vel;
    if(i1 > i2) boundary_vel = std::max(0.0, vel_x + boost_vel);
    else { assert(i2 > i1); boundary_vel = std::min(0.0, vel_x + boost_vel); }
    double ghost_vel = (dist*(boundary_vel) - 0.5*grids[d_x](i2,j2)*vel_x)/(0.5*grids[d_x](i3,j3));
    m_mom_x(i1,j1) = m_rho(i1,j1)*ghost_vel;
    m_mom_x(i2,j2) = m_rho(i2,j2)*ghost_vel;
    m_mom_y(i1,j1) = m_rho(i1,j1)*vel_y;
    m_mom_y(i2,j2) = m_rho(i2,j2)*vel_y;
  }
}

//Applies the wall boundary condition to the cells indexed by the given indices
//assuming that (i1,j1) is at the edge of the domain and (i4,j4) is on the interior
//Assumes two ghost zones (i.e. N_GHOST == 2), will abort if not
//Meant to be called repeatedly during advanceTime; energy is handled, not temperature
void PlasmaDomain::reflectBoundaryExtrapolate(std::vector<Grid> &grids, int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  assert(N_GHOST == 2 && "This function assumes two ghost zones");
  Grid &m_rho = grids[rho], &m_thermal_energy = grids[thermal_energy], &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y];
  m_rho(i1,j1) = m_rho(i3,j3); m_rho(i2,j2) = m_rho(i3,j3); //Match density of nearest interior cell
  m_thermal_energy(i1,j1) = m_thermal_energy(i3,j3); m_thermal_energy(i2,j2) = m_thermal_energy(i3,j3); //Match temperature of nearest interior cell
  m_mom_x(i1,j1) = 0.0; m_mom_x(i2,j2) = 0.0; m_mom_x(i3,j3) = 0.0; //Halt all advection at the wall boundary
  m_mom_y(i1,j1) = 0.0; m_mom_y(i2,j2) = 0.0; m_mom_y(i3,j3) = 0.0;
}

//Applies the wall boundary condition to the cells indexed by the given indices
//assuming that (i1,j1) is at the edge of the domain and (i4,j4) is on the interior
//Assumes two ghost zones (i.e. N_GHOST == 2), will abort if not
//Meant to be called repeatedly during advanceTime; energy is handled, not temperature
void PlasmaDomain::fixedBoundaryExtrapolate(std::vector<Grid> &grids, int i1, int i2, int i3, int i4, int j1, int j2, int j3, int j4)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  // Currently does nothing; leaves the ghost cells at their initial values
  // assert(N_GHOST == 2 && "This function assumes two ghost zones");
  Grid &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y];
  m_mom_x(i1,j1) = 0.0; m_mom_x(i2,j2) = 0.0; m_mom_x(i3,j3) = 0.0; //Halt all advection at the wall boundary
  m_mom_y(i1,j1) = 0.0; m_mom_y(i2,j2) = 0.0; m_mom_y(i3,j3) = 0.0;
}

// Applies two-dimensional Sovitzky-Golay filter to member variable m_grids
// with 3x3 window and 2x2-order fitting polynomials
void PlasmaDomain::filterSavitzkyGolay()
{
  filterSavitzkyGolay(m_grids);
}

// Applies two-dimensional Sovitzky-Golay filter to given Grid vector
// with 5x5 window and 3x3-order fitting polynomials
void PlasmaDomain::filterSavitzkyGolay(std::vector<Grid> &grids)
{
  assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
  std::vector<Grid> filtered = grids;
  #pragma omp parallel
  {
    #pragma omp for
    for(int varname : {rho, thermal_energy}){
      singleVarSavitzkyGolay(grids[varname]);
    }
  }
  propagateChanges(grids);
}

void PlasmaDomain::singleVarSavitzkyGolay(Grid &grid)
{
  assert(N_GHOST == 2 && "S-G filtering implementation assumes two ghost zones");
  static const double coeff[] =
    {+7.346939E-03,	-2.938776E-02,	-4.163265E-02,	-2.938776E-02,	+7.346939E-03,
    -2.938776E-02,	+1.175510E-01,	+1.665306E-01,	+1.175510E-01,	-2.938776E-02,
    -4.163265E-02,	+1.665306E-01,	+2.359184E-01,  +1.665306E-01,  -4.163265E-02,
    -2.938776E-02,  +1.175510E-01,  +1.665306E-01,  +1.175510E-01,  -2.938776E-02,
    +7.346939E-03,  -2.938776E-02,  -4.163265E-02,  -2.938776E-02,  +7.346939E-03}; //from Chandra Sekhar, 2015
  int xdim = grid.rows();
  int ydim = grid.cols();
  Grid filtered = grid;
  for (int center_i = m_xl; center_i <= m_xu; center_i++){
    for(int center_j = m_yl; center_j <= m_yu; center_j++){
      int i[] = {center_i-2, center_i-1, center_i, center_i+1, center_i+2};
      int j[] = {center_j-2, center_j-1, center_j, center_j+1, center_j+2};
      if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
        i[0] = (i[0]+xdim)%xdim;
        i[1] = (i[1]+xdim)%xdim;
        i[3] = (i[3]+xdim)%xdim;
        i[4] = (i[4]+xdim)%xdim;
      }
      if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
        j[0] = (j[0]+ydim)%ydim;
        j[1] = (j[1]+ydim)%ydim;
        j[3] = (j[3]+ydim)%ydim;
        j[4] = (j[4]+ydim)%ydim;
      }
      double filtered_val = 0.0;
      // double total_weight = 0.0;
      for(int u = 0; u < 5; u++){
        for(int v = 0; v < 5; v++){
          int curr_i = i[u], curr_j = j[v];
          // if(curr_i >= m_xl && curr_i <= m_xu && curr_j >= m_yl && curr_j <= m_yu){
            filtered_val = filtered_val + coeff[u*5 + v]*grid(curr_j,curr_j);
            // total_weight = total_weight + coeff[u*5 + v];
          // }
        }
      }
      // filtered_val = filtered_val/total_weight;
      filtered(center_i,center_j) = filtered_val;
    }
  }
  grid = filtered;
}