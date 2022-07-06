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
    Grid i_viscous_force_x = visc_coeff*m_pd.laplacian(grids[i_mom_x]);
    Grid i_viscous_force_y = visc_coeff*m_pd.laplacian(grids[i_mom_y]);
    Grid e_viscous_force_x = visc_coeff*m_pd.laplacian(grids[e_mom_x]);
    Grid e_viscous_force_y = visc_coeff*m_pd.laplacian(grids[e_mom_y]);
    // lorentz forces
    std::vector<Grid> b = {grids[b_x],grids[b_y]};
    Grid i_v_cross_B = Grid::CrossProduct2D(v_i,b);
    Grid e_v_cross_B = Grid::CrossProduct2D(v_e,b);
    Grid i_F_x = E*grids[i_n]*(grids[E_x]+i_v_cross_B);
    Grid i_F_y = E*grids[i_n]*(grids[E_y]+i_v_cross_B);
    Grid e_F_x = -E*grids[e_n]*(grids[E_x]+e_v_cross_B);
    Grid e_F_y = -E*grids[e_n]*(grids[E_y]+e_v_cross_B); 
    // momentum equations   
    Grid i_d_mom_x_dt = -m_pd.transportDivergence2D(grids[i_mom_x], v_i)
                        - m_pd.derivative1D(grids[i_press], 0)
                        + m_pd.m_ghost_zone_mask * (grids[i_rho]*grids[grav_x] + i_viscous_force_x + i_F_x);
    Grid i_d_mom_y_dt = -m_pd.transportDivergence2D(grids[i_mom_y], v_i)
                        - m_pd.derivative1D(grids[i_press], 1)
                        + m_pd.m_ghost_zone_mask * (grids[i_rho]*grids[grav_y] + i_viscous_force_y + i_F_y);
    Grid e_d_mom_x_dt = -m_pd.transportDivergence2D(grids[e_mom_x], v_e)
                        - m_pd.derivative1D(grids[e_press], 0)
                        + m_pd.m_ghost_zone_mask * (grids[e_rho]*grids[grav_x] + e_viscous_force_x + e_F_x);
    Grid e_d_mom_y_dt = -m_pd.transportDivergence2D(grids[e_mom_y], v_e)
                        - m_pd.derivative1D(grids[e_press], 1)
                        + m_pd.m_ghost_zone_mask * (grids[e_rho]*grids[grav_y] + e_viscous_force_y + e_F_y);
    // energy equations
    Grid i_d_thermal_dt =  - m_pd.transportDivergence2D(grids[i_thermal_energy],v_i)
                                    - grids[i_press]*m_pd.divergence2D(v_i);
    Grid e_d_thermal_dt =  - m_pd.transportDivergence2D(grids[e_thermal_energy],v_e)
                                    - grids[e_press]*m_pd.divergence2D(v_e);
    // induction equations
    std::vector<Grid> curl_E = m_pd.curlZ(grids[E_z]);
    Grid d_bi_x_dt = -curl_E[0];
    Grid d_bi_y_dt = -curl_E[1];
    // return time derivatives
    return {i_d_rho_dt,e_d_rho_dt,i_d_mom_x_dt,i_d_mom_y_dt,e_d_mom_x_dt,e_d_mom_y_dt,i_d_thermal_dt,e_d_thermal_dt,d_bi_x_dt,d_bi_y_dt};
}

void Ideal2F::computeConstantTerms(std::vector<Grid> &grids){
}

void Ideal2F::recomputeDerivedVariables(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    recomputeTemperature(grids);
    recomputeKineticEnergy(grids);
    recomputeMagneticFields(grids);
    recomputeElectricFields(grids);
    recomputeDT();
}

void Ideal2F::recomputeTemperature(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    // enforce minimum on thermal energies
    grids[i_thermal_energy] = grids[i_thermal_energy].max(m_pd.thermal_energy_min);
    grids[e_thermal_energy] = grids[e_thermal_energy].max(m_pd.thermal_energy_min);
    // recompute pressures
    grids[i_press] = (m_pd.m_adiabatic_index - 1.0)*grids[i_thermal_energy];
    grids[e_press] = (m_pd.m_adiabatic_index - 1.0)*grids[e_thermal_energy];
    grids[press] = grids[i_press] + grids[e_press];
    // recompute temperatures and enforce minimum
    grids[i_temp] = (grids[i_press]/(K_B*grids[i_n])).max(m_pd.temp_min);
    grids[e_temp] = (grids[e_press]/(K_B*grids[e_n])).max(m_pd.temp_min);
}

void Ideal2F::recomputeThermalEnergy(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    grids[i_press] = K_B*grids[i_n]*grids[i_temp];
    grids[e_press] = K_B*grids[e_n]*grids[e_temp];
    grids[press] = grids[i_press] + grids[e_press];
    grids[i_thermal_energy] = grids[i_press]/(m_pd.m_adiabatic_index - 1.0);
    grids[e_thermal_energy] = grids[e_press]/(m_pd.m_adiabatic_index - 1.0);
}

void Ideal2F::recomputeKineticEnergy(std::vector<Grid> &grids){
    grids[i_v_x] = grids[i_mom_x]/grids[i_rho];
    grids[i_v_y] = grids[i_mom_y]/grids[i_rho];
    grids[e_v_x] = grids[e_mom_x]/grids[e_rho];
    grids[e_v_y] = grids[e_mom_y]/grids[e_rho];
    grids[i_kinetic_energy] = 0.5*grids[i_rho]*(grids[i_v_x].square()+grids[i_v_y].square());
    grids[e_kinetic_energy] = 0.5*grids[e_rho]*(grids[e_v_x].square()+grids[e_v_y].square());
}

void Ideal2F::recomputeNumberDensity(std::vector<Grid> &grids){
    grids[i_n] = grids[i_rho]/m_pd.m_ion_mass;
    grids[e_n] = grids[e_rho]/M_ELECTRON;
}

void Ideal2F::recomputeMagneticFields(std::vector<Grid> &grids){
    grids[b_x] = m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x];
    grids[b_y] = m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_mag] = (grids[b_x].square() + grids[b_y].square()).sqrt();
    grids[b_hat_x] = grids[b_x]/grids[b_mag];
    grids[b_hat_y] = grids[b_y]/grids[b_mag];
    // ensure bhat is zero when b is zero
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(grids[b_mag](i,j) == 0.0){
                grids[b_hat_x](i,j) = 0.0;
                grids[b_hat_y](i,j) = 0.0;
            }
        }
    }
}

void Ideal2F::recomputeElectricFields(std::vector<Grid> &grids){
    grids[E_x] = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    grids[E_y] = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    grids[E_z] = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
}

void Ideal2F::catchUnderdensity(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    for (int ind : densities()){
        for(int i=0; i<m_pd.m_xdim; i++){
            for(int j=0; j<m_pd.m_ydim; j++){
                if(grids[ind](i,j) < m_pd.rho_min){
                    if(i >= m_pd.m_xl && i <= m_pd.m_xu && j >= m_pd.m_yl && j <= m_pd.m_yu){
                        grids[ind](i,j) = m_pd.rho_min;
                    }
                }
            }
        }
    }
}

void Ideal2F::recomputeDT(){
    // sound speed for a neutral plasma - intertia carried by ions
    Grid c_s = (m_pd.m_adiabatic_index*m_grids[press]/m_grids[i_rho]).sqrt();
    Grid c_s_sq = c_s.square();
    // alfen speed for ion component
    Grid v_alfven = m_grids[b_mag]/(4.0*PI*m_grids[i_rho]).sqrt();
    Grid v_alfven_sq = v_alfven.square();
    // magnetoacoustic modes
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);
    Grid delta = (one - 4.0*c_s_sq*v_alfven_sq/(c_s_sq+v_alfven_sq).square()).sqrt();
    Grid v_fast = (0.5*(c_s_sq + v_alfven_sq)*(one + delta)).sqrt();
    Grid v_slow = (0.5*(c_s_sq + v_alfven_sq)*(one - delta)).sqrt();
    // velocity magnitude
    Grid vel_mag = (m_grids[e_v_x].square() + m_grids[e_v_y].square()).sqrt();
    // get timestep
    Grid diagonals = (m_pd.m_grids[PlasmaDomain::d_x].square() + m_pd.m_grids[PlasmaDomain::d_y].square()).sqrt();
    m_grids[dt] = diagonals/(c_s + v_alfven + v_fast + v_slow + vel_mag);
}

Grid Ideal2F::getDT(){
    return m_grids[dt];
}

void Ideal2F::propagateChanges(std::vector<Grid> &grids)
{
    catchUnderdensity(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVariables(grids);
}

void Ideal2F::populateVariablesFromState(std::vector<Grid> &grids){
    computeConstantTerms(grids);
    recomputeThermalEnergy(grids);
    recomputeDerivedVariables(grids);
    m_pd.updateGhostZones();
}