#include "idealmhdcons.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>

IdealMHDCons::IdealMHDCons(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void IdealMHDCons::parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){}

void IdealMHDCons::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[rho] += step*(m_pd.m_ghost_zone_mask*time_derivatives[0]);
    grids[mom_x] += step*(m_pd.m_ghost_zone_mask*time_derivatives[1]);
    grids[mom_y] += step*(m_pd.m_ghost_zone_mask*time_derivatives[2]);
    grids[energy] += step*(m_pd.m_ghost_zone_mask*time_derivatives[3]);
    grids[bi_x] += step*(m_pd.m_ghost_zone_mask*time_derivatives[4]);
    grids[bi_y] += step*(m_pd.m_ghost_zone_mask*time_derivatives[5]);
    propagateChanges(grids);
}

std::vector<Grid> IdealMHDCons::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // PlasmaDomain grid references for more concise notation
    Grid& be_x = m_pd.m_grids[PlasmaDomain::be_x];
    Grid& be_y = m_pd.m_grids[PlasmaDomain::be_y];
    // continuity equations
    std::vector<Grid> v = {grids[v_x],grids[v_y]};
    Grid d_rho_dt = -m_pd.transportDivergence2D(grids[rho],v);
    // viscous forces
    Grid viscous_force_x = visc_coeff*grids[rho]*m_pd.laplacian(grids[v_x]);
    Grid viscous_force_y = visc_coeff*grids[rho]*m_pd.laplacian(grids[v_y]);
    // tensor for momentum equation
    Grid Txx = grids[press_tot] + grids[b_x]*grids[b_x]/(4.*PI);
    Grid Tyy = grids[press_tot] + grids[b_y]*grids[b_y]/(4.*PI);
    Grid Txy = grids[b_x]*grids[b_y]/(4.*PI);
    Grid Tyx = grids[b_y]*grids[b_x]/(4.*PI);
    // momentum equations   
    Grid d_mom_x_dt =   - m_pd.transportDivergence2D(grids[mom_x], v)
                        + m_pd.m_ghost_zone_mask * (grids[rho]*grids[grav_x] + viscous_force_x - m_pd.divergence2D({Txx,Tyx}));
    Grid d_mom_y_dt =   - m_pd.transportDivergence2D(grids[mom_y], v)
                        + m_pd.m_ghost_zone_mask * (grids[rho]*grids[grav_y] + viscous_force_y - m_pd.divergence2D({Txy,Tyy}));
    // energy equations
    Grid v_dot_B = Grid::DotProduct2D(v,{grids[b_x],grids[b_y]});
    Grid Ux = v_dot_B*grids[b_x]/(4.*PI)-grids[press_tot]*grids[v_x];
    Grid Uy = v_dot_B*grids[b_y]/(4.*PI)-grids[press_tot]*grids[v_y];
    Grid d_energy_dt =  - m_pd.transportDivergence2D(grids[energy],v)
                        + m_pd.divergence2D({Ux,Uy});
    // induction equations
    Grid Yxx = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid Yyy = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid Yxy = grids[v_x]*grids[b_y] - grids[b_x]*grids[v_y];
    Grid Yyx = -1.*Yxy;
    Grid d_bi_x_dt = -m_pd.divergence2D({Yxx,Yyx});
    Grid d_bi_y_dt = -m_pd.divergence2D({Yyy,Yxy});
    // return time derivatives
    return {d_rho_dt, d_mom_x_dt, d_mom_y_dt, d_energy_dt, d_bi_x_dt, d_bi_y_dt};
}

Grid IdealMHDCons::getDT(){
    return m_grids[dt];
}

void IdealMHDCons::populateVariablesFromState(std::vector<Grid> &grids){
    computeConstantTerms(grids);
    recomputeKineticEnergy(grids);
    recomputeMagneticEnergy(grids);
    recomputeThermalEnergyFromTemp(grids);
    grids[energy] = grids[thermal_energy] + grids[kinetic_energy] + grids[mag_energy];
    recomputeDerivedVariables(grids);
    m_pd.updateGhostZones();
}

void IdealMHDCons::computeConstantTerms(std::vector<Grid> &grids){
}

void IdealMHDCons::propagateChanges(std::vector<Grid> &grids)
{
    enforceMinimums(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVariables(grids);
}

void IdealMHDCons::recomputeDerivedVariables(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    recomputeKineticEnergy(grids);
    recomputeMagneticEnergy(grids);
    recomputeThermalEnergy(grids);
    recomputeTemperature(grids);
    recomputeDT();
}

void IdealMHDCons::recomputeNumberDensity(std::vector<Grid> &grids){
    grids[n] = grids[rho]/m_pd.m_ion_mass;
}

void IdealMHDCons::recomputeKineticEnergy(std::vector<Grid> &grids){
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[kinetic_energy] = 0.5*grids[rho]*(grids[v_x].square()+grids[v_y].square());
}

void IdealMHDCons::recomputeMagneticEnergy(std::vector<Grid> &grids){
    grids[b_x] = m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x];
    grids[b_y] = m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_mag] = (grids[b_x].square() + grids[b_y].square()).sqrt();
    grids[mag_energy] = grids[b_mag].square()/(8.*PI);
    // compute field direction and ensure is zero when B=0
    grids[b_hat_x] = grids[b_x]/grids[b_mag];
    grids[b_hat_y] = grids[b_y]/grids[b_mag];
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(grids[b_mag](i,j) == 0.0){
                grids[b_hat_x](i,j) = 0.0;
                grids[b_hat_y](i,j) = 0.0;
            }
        }
    }
}

void IdealMHDCons::recomputeThermalEnergy(std::vector<Grid> &grids)
{
    grids[thermal_energy] = grids[energy] - grids[kinetic_energy] - grids[mag_energy];
}

void IdealMHDCons::recomputeTemperature(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    // enforce minimum on thermal energies
    grids[thermal_energy] = grids[thermal_energy].max(m_pd.thermal_energy_min);
    // recompute pressures
    grids[press] = (m_pd.m_adiabatic_index - 1.0)*grids[thermal_energy];
    grids[press_tot] = grids[press] + grids[mag_energy];
    // recompute temperatures and enforce minimum
    grids[temp] = (grids[press]/(2*K_B*grids[n])).max(m_pd.temp_min);
}

void IdealMHDCons::recomputeThermalEnergyFromTemp(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    grids[press] = 2*K_B*grids[n]*grids[temp];
    grids[thermal_energy] = grids[press]/(m_pd.m_adiabatic_index - 1.0);
}

void IdealMHDCons::recomputeDT(){
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

void IdealMHDCons::enforceMinimums(std::vector<Grid>& grids)
{
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[rho] = grids[n]*m_pd.m_ion_mass;
    grids[thermal_energy] = grids[thermal_energy].max(m_pd.thermal_energy_min);
}