#include "ideal2F.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>

Ideal2F::Ideal2F(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void Ideal2F::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[i_rho] += step*m_pd.m_ghost_zone_mask*time_derivatives[0];
    grids[e_rho] += step*m_pd.m_ghost_zone_mask*time_derivatives[1];
    grids[i_mom_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[2];
    grids[i_mom_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[3];
    grids[e_mom_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[4];
    grids[e_mom_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[5];
    grids[i_thermal_energy] += step*m_pd.m_ghost_zone_mask*time_derivatives[6];
    grids[e_thermal_energy] += step*m_pd.m_ghost_zone_mask*time_derivatives[7];
    grids[bi_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[8];
    grids[bi_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[9];
    grids[bi_z] += step*m_pd.m_ghost_zone_mask*time_derivatives[10];
    grids[E_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[11];
    grids[E_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[12];
    grids[E_z] += step*m_pd.m_ghost_zone_mask*time_derivatives[13];
    propagateChanges(grids);
}

std::vector<Grid> Ideal2F::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // continuity equations
    std::vector<Grid> v_i = {grids[i_v_x],grids[i_v_y]};
    std::vector<Grid> v_e = {grids[e_v_x],grids[e_v_y]};
    Grid i_d_rho_dt = -m_pd.transportDivergence2D(grids[i_rho],v_i);
    Grid e_d_rho_dt = -m_pd.transportDivergence2D(grids[e_rho],v_e);
    // viscous forces
    Grid i_viscous_force_x = visc_coeff*grids[i_rho]*m_pd.laplacian(grids[i_v_x]);
    Grid i_viscous_force_y = visc_coeff*grids[i_rho]*m_pd.laplacian(grids[i_v_y]);
    Grid e_viscous_force_x = visc_coeff*grids[e_rho]*m_pd.laplacian(grids[e_v_x]);
    Grid e_viscous_force_y = visc_coeff*grids[e_rho]*m_pd.laplacian(grids[e_v_y]);
    // lorentz forces
    std::vector<Grid> b = {grids[b_x],grids[b_y]};
    std::vector<Grid> i_v_cross_B(2,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    i_v_cross_B[0] = grids[i_v_y]*grids[b_z];
    i_v_cross_B[1] = -grids[i_v_x]*grids[b_z];
    std::vector<Grid> e_v_cross_B(2,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    e_v_cross_B[0] = grids[e_v_y]*grids[b_z];
    e_v_cross_B[1] = -grids[e_v_x]*grids[b_z];
    Grid i_F_x =  E*grids[i_n]*(grids[E_x] + i_v_cross_B[0]);
    Grid i_F_y =  E*grids[i_n]*(grids[E_y] + i_v_cross_B[1]);
    Grid e_F_x = -E*grids[e_n]*(grids[E_x] + e_v_cross_B[0]);
    Grid e_F_y = -E*grids[e_n]*(grids[E_y] + e_v_cross_B[1]); 
    // momentum equations   
    Grid i_d_mom_x_dt = - m_pd.transportDivergence2D(grids[i_mom_x], v_i)
                        - m_pd.derivative1D(grids[i_press], 0)
                        + m_pd.m_ghost_zone_mask * (grids[i_rho]*grids[grav_x] + i_viscous_force_x + i_F_x);
    Grid i_d_mom_y_dt = - m_pd.transportDivergence2D(grids[i_mom_y], v_i)
                        - m_pd.derivative1D(grids[i_press], 1)
                        + m_pd.m_ghost_zone_mask * (grids[i_rho]*grids[grav_y] + i_viscous_force_y + i_F_y);
    Grid e_d_mom_x_dt = - m_pd.transportDivergence2D(grids[e_mom_x], v_e)
                        - m_pd.derivative1D(grids[e_press], 0)
                        + m_pd.m_ghost_zone_mask * (grids[e_rho]*grids[grav_x] + e_viscous_force_x + e_F_x);
    Grid e_d_mom_y_dt = - m_pd.transportDivergence2D(grids[e_mom_y], v_e)
                        - m_pd.derivative1D(grids[e_press], 1)
                        + m_pd.m_ghost_zone_mask * (grids[e_rho]*grids[grav_y] + e_viscous_force_y + e_F_y);
    // energy equations
    Grid i_d_thermal_dt =   - m_pd.transportDivergence2D(grids[i_thermal_energy],v_i)
                            - grids[i_press]*m_pd.divergence2D(v_i);
    Grid e_d_thermal_dt =   - m_pd.transportDivergence2D(grids[e_thermal_energy],v_e)
                            - grids[e_press]*m_pd.divergence2D(v_e);
    // magnetic field propagation
    Grid d_bi_x_dt = -m_pd.derivative1D(grids[E_z],1);
    Grid d_bi_y_dt =  m_pd.derivative1D(grids[E_z],0);
    Grid d_bi_z_dt = m_pd.derivative1D(grids[E_x],1) - m_pd.derivative1D(grids[E_y],0);
    // electric field propagation
    Grid d_E_x_dt =  C*m_pd.derivative1D(grids[b_z],1) - 4*PI*grids[j_x];
    Grid d_E_y_dt = -C*m_pd.derivative1D(grids[b_z],0) - 4*PI*grids[j_y];
    Grid d_E_z_dt =  C*(m_pd.derivative1D(grids[b_y],0) - m_pd.derivative1D(grids[b_x],1));
    // return time derivatives
    return {i_d_rho_dt,e_d_rho_dt,i_d_mom_x_dt,i_d_mom_y_dt,e_d_mom_x_dt,e_d_mom_y_dt,i_d_thermal_dt,e_d_thermal_dt,
        d_bi_x_dt,d_bi_y_dt,d_bi_z_dt,d_E_x_dt,d_E_y_dt,d_E_z_dt};
}

void Ideal2F::populateVariablesFromState(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeEvolvedVarsFromStateVars(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVarsFromEvolvedVars(grids);
    recomputeDT();
}

void Ideal2F::propagateChanges(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    enforceMinimums(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVarsFromEvolvedVars(grids);
    recomputeDT();
}

void Ideal2F::enforceMinimums(std::vector<Grid>& grids)
{
    grids[i_n] = (grids[i_rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[e_n] = (grids[e_rho]/M_ELECTRON).max(m_pd.density_min);
    grids[i_rho] = grids[i_n]*m_pd.m_ion_mass;
    grids[e_rho] = grids[e_n]*M_ELECTRON;
    grids[i_thermal_energy] = grids[i_thermal_energy].max(m_pd.thermal_energy_min);
    grids[e_thermal_energy] = grids[e_thermal_energy].max(m_pd.thermal_energy_min);
}

void Ideal2F::recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[i_n] = grids[i_rho]/m_pd.m_ion_mass;
    grids[e_n] = grids[e_rho]/M_ELECTRON;
    
    grids[i_press] = grids[i_n]*K_B*grids[i_temp];
    grids[e_press] = grids[e_n]*K_B*grids[e_temp];
    grids[press] = grids[i_press] + grids[e_press];
    grids[i_thermal_energy] = grids[i_press]/(m_pd.m_adiabatic_index - 1.0);
    grids[e_thermal_energy] = grids[e_press]/(m_pd.m_adiabatic_index - 1.0);
    grids[E_x] = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    grids[E_y] = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    grids[E_z] = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
}

void Ideal2F::recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[i_n] = grids[i_rho]/m_pd.m_ion_mass;
    grids[e_n] = grids[e_rho]/M_ELECTRON;
    grids[i_v_x] = grids[i_mom_x]/grids[i_rho];
    grids[i_v_y] = grids[i_mom_y]/grids[i_rho];
    grids[e_v_x] = grids[e_mom_x]/grids[e_rho];
    grids[e_v_y] = grids[e_mom_y]/grids[e_rho];
    grids[j_x] = grids[i_n]*E*grids[i_v_x] - grids[e_n]*E*grids[e_v_x];
    grids[j_y] = grids[i_n]*E*grids[i_v_y] - grids[e_n]*E*grids[e_v_y];
    grids[i_press] = (m_pd.m_adiabatic_index - 1.0)*grids[i_thermal_energy];
    grids[e_press] = (m_pd.m_adiabatic_index - 1.0)*grids[e_thermal_energy];
    grids[press] = grids[i_press] + grids[e_press];
    grids[i_temp] = (grids[i_press]/(K_B*grids[i_n])).max(m_pd.temp_min);
    grids[e_temp] = (grids[e_press]/(K_B*grids[e_n])).max(m_pd.temp_min);
    grids[b_x] = m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x];
    grids[b_y] = m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_z] = grids[bi_z]; // code assumes be_z = 0 for now
    grids[b_mag] = (grids[b_x].square() + grids[b_y].square() + grids[b_z].square()).sqrt();
    grids[b_mag_xy] = (grids[b_x].square() + grids[b_y].square()).sqrt();
    grids[b_hat_x] = grids[b_x]/grids[b_mag_xy];
    grids[b_hat_y] = grids[b_y]/grids[b_mag_xy];
    catchNullFieldDirection(grids);
    grids[rho] = grids[i_rho] + grids[e_rho];
    grids[rho_c] = E*(grids[i_n] - grids[e_n]);
    grids[n] = grids[i_rho]/m_pd.m_ion_mass;
    grids[dn] = grids[i_n] - grids[e_n];
    grids[divBcond] = m_pd.divergence2D(grids[b_x],grids[b_y]);
    grids[divEcond] = m_pd.divergence2D(grids[E_x],grids[E_y]) - 4*PI*grids[rho_c];
}

void Ideal2F::catchNullFieldDirection(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(grids[b_mag_xy](i,j) == 0.0){
                grids[b_hat_x](i,j) = 0.0;
                grids[b_hat_y](i,j) = 0.0;
            }
        }
    }
}

void Ideal2F::recomputeDT(){
    Grid vel_mag = (m_grids[e_v_x].square() + m_grids[e_v_y].square()).sqrt();
    Grid vel_th = (K_B*m_grids[e_temp]/M_ELECTRON).sqrt();
    Grid w_pe = (4*PI*m_grids[e_n]*E*E/M_ELECTRON).sqrt();
    Grid tau_pe = 1/w_pe;
    Grid diagonals = (m_pd.m_grids[PlasmaDomain::d_x].square() + m_pd.m_grids[PlasmaDomain::d_y].square()).sqrt();
    m_grids[dt] = diagonals/(3*vel_th + vel_mag);
    if (tau_pe.min() < m_grids[dt].min()) m_grids[dt] = tau_pe;
}