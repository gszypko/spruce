#include "idealmhd2E.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>

IdealMHD2E::IdealMHD2E(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void IdealMHD2E::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[rho] += step*m_pd.m_ghost_zone_mask*time_derivatives[0];
    grids[mom_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[1];
    grids[mom_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[2];
    grids[i_thermal_energy] += step*m_pd.m_ghost_zone_mask*time_derivatives[3];
    grids[e_thermal_energy] += step*m_pd.m_ghost_zone_mask*time_derivatives[4];
    grids[bi_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[5];
    grids[bi_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[6];
    propagateChanges(grids);
}

std::vector<Grid> IdealMHD2E::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // PlasmaDomain grid references for more concise notation
    Grid& be_x = m_pd.m_grids[PlasmaDomain::be_x];
    Grid& be_y = m_pd.m_grids[PlasmaDomain::be_y];
    // continuity equations
    std::vector<Grid> v = {grids[v_x],grids[v_y]};
    Grid d_rho_dt = -m_pd.transportDivergence2D(grids[rho],v);
    // viscous forces
    Grid viscous_force_x = visc_coeff*m_pd.laplacian(grids[mom_x]);
    Grid viscous_force_y = visc_coeff*m_pd.laplacian(grids[mom_y]);
    // magnetic forces
    Grid curl_db = m_pd.curl2D(grids[bi_x],grids[bi_y])/(4.0*PI);
    std::vector<Grid> external_mag_force = Grid::CrossProductZ2D(curl_db,{be_x,be_y});
    std::vector<Grid> internal_mag_force = Grid::CrossProductZ2D(curl_db,{grids[bi_x],grids[bi_y]}); //second order terms
    // momentum equations   
    Grid d_mom_x_dt =   - m_pd.transportDivergence2D(grids[mom_x], v)
                        - m_pd.derivative1D(grids[press], 0)
                        + m_pd.m_ghost_zone_mask * (grids[rho]*grids[grav_x] + viscous_force_x
                        + external_mag_force[0] + internal_mag_force[0]);
    Grid d_mom_y_dt =   - m_pd.transportDivergence2D(grids[mom_y], v)
                        - m_pd.derivative1D(grids[press], 1)
                        + m_pd.m_ghost_zone_mask * (grids[rho]*grids[grav_y] + viscous_force_y
                        + external_mag_force[1] + internal_mag_force[1]);
    // energy equations
    Grid i_d_thermal_energy_dt =   - m_pd.transportDivergence2D(grids[i_thermal_energy],v)
                            - grids[i_press]*m_pd.divergence2D(v);
    Grid e_d_thermal_energy_dt =   - m_pd.transportDivergence2D(grids[e_thermal_energy],v)
                            - grids[e_press]*m_pd.divergence2D(v);
    // induction equations
    std::vector<Grid> induction_rhs_external = m_pd.curlZ(Grid::CrossProduct2D(v,{be_x,be_y}));
    std::vector<Grid> induction_rhs_internal = m_pd.curlZ(Grid::CrossProduct2D(v,{grids[bi_x],grids[bi_y]}));
    Grid d_bi_x_dt = induction_rhs_external[0] + induction_rhs_internal[0];
    Grid d_bi_y_dt = induction_rhs_external[1] + induction_rhs_internal[1];
    // return time derivatives
    return {d_rho_dt, d_mom_x_dt, d_mom_y_dt, i_d_thermal_energy_dt, e_d_thermal_energy_dt, d_bi_x_dt, d_bi_y_dt};
}

void IdealMHD2E::computeConstantTerms(std::vector<Grid> &grids){
}

void IdealMHD2E::recomputeDerivedVariables(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    recomputeTemperature(grids);
    recomputeKineticEnergy(grids);
    recomputeMagneticFields(grids);
    recomputeDT();
}

void IdealMHD2E::recomputeTemperature(std::vector<Grid> &grids){
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
    grids[i_temp] = (grids[i_press]/(K_B*grids[n])).max(m_pd.temp_min);
    grids[e_temp] = (grids[e_press]/(K_B*grids[n])).max(m_pd.temp_min);
}

void IdealMHD2E::recomputeThermalEnergy(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    grids[i_press] = K_B*grids[n]*grids[i_temp];
    grids[e_press] = K_B*grids[n]*grids[e_temp];
    grids[press] = grids[i_press] + grids[e_press];
    grids[i_thermal_energy] = grids[i_press]/(m_pd.m_adiabatic_index - 1.0);
    grids[e_thermal_energy] = grids[e_press]/(m_pd.m_adiabatic_index - 1.0);
}

void IdealMHD2E::recomputeKineticEnergy(std::vector<Grid> &grids){
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[kinetic_energy] = 0.5*grids[rho]*(grids[v_x].square()+grids[v_y].square());
}

void IdealMHD2E::recomputeNumberDensity(std::vector<Grid> &grids){
    grids[n] = grids[rho]/m_pd.m_ion_mass;
}

void IdealMHD2E::recomputeMagneticFields(std::vector<Grid> &grids){
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

void IdealMHD2E::recomputeDT(){
    Grid c_s = (m_pd.m_adiabatic_index*m_grids[press]/m_grids[rho]).sqrt();
    Grid c_s_sq = c_s.square();

    Grid v_alfven = m_grids[b_mag]/(4.0*PI*m_grids[rho]).sqrt();
    Grid v_alfven_sq = v_alfven.square();
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);

    //Magnetoacoustic modes
    Grid delta = (one - 4.0*c_s_sq*v_alfven_sq/(c_s_sq+v_alfven_sq).square()).sqrt();
    Grid v_fast = (0.5*(c_s_sq + v_alfven_sq)*(one + delta)).sqrt();
    Grid v_slow = (0.5*(c_s_sq + v_alfven_sq)*(one - delta)).sqrt();

    //Bulk velocity transit time
    Grid diagonals = (m_pd.m_grids[PlasmaDomain::d_x].square() + m_pd.m_grids[PlasmaDomain::d_y].square()).sqrt();
    Grid vel_mag = (m_grids[v_x].square() + m_grids[v_y].square()).sqrt();

    m_grids[dt] = diagonals/(c_s + v_alfven + v_fast + v_slow + vel_mag);
}

Grid IdealMHD2E::getDT(){
    return m_grids[dt];
}

void IdealMHD2E::propagateChanges(std::vector<Grid> &grids)
{
    enforceMinimums(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVariables(grids);
}

void IdealMHD2E::populateVariablesFromState(std::vector<Grid> &grids){
    computeConstantTerms(grids);
    recomputeThermalEnergy(grids);
    recomputeDerivedVariables(grids);
    m_pd.updateGhostZones();
}

void IdealMHD2E::enforceMinimums(std::vector<Grid>& grids)
{
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[rho] = grids[n]*m_pd.m_ion_mass;
    grids[i_thermal_energy] = grids[i_thermal_energy].max(m_pd.thermal_energy_min);
    grids[e_thermal_energy] = grids[e_thermal_energy].max(m_pd.thermal_energy_min);
}