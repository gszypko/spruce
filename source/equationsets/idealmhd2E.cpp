#include "idealmhd2E.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>

IdealMHD2E::IdealMHD2E(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void IdealMHD2E::parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for (int i=0; i<lhs.size(); i++){
        if (lhs[i] == "global_viscosity") m_global_viscosity = stod(rhs[i]);
        else if (lhs[i] == "viscosity_opt") m_viscosity_opt = rhs[i];
        else{
            std::cerr << lhs[i] << " is not recognized for this equation set." << std::endl;
            assert(false);
        }
    }
}

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

std::vector<Grid> IdealMHD2E::computeTimeDerivatives(const std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // PlasmaDomain grid references for more concise notation
    Grid& d_x = m_pd.m_grids[PlasmaDomain::d_x];
    Grid& d_y = m_pd.m_grids[PlasmaDomain::d_y];
    Grid& be_x = m_pd.m_grids[PlasmaDomain::be_x];
    Grid& be_y = m_pd.m_grids[PlasmaDomain::be_y];
    // continuity equations
    std::vector<Grid> v = {grids[v_x],grids[v_y]};
    Grid d_rho_dt = -m_pd.transportDivergence2D(grids[rho],v);
    // viscous forces
    Grid global_visc_coeff = m_global_viscosity*0.5*(d_x.square()+d_y.square())/grids[dt];
    Grid viscous_force_x, viscous_force_y;
    if (m_viscosity_opt == "momentum"){
        viscous_force_x = global_visc_coeff*m_pd.laplacian(grids[mom_x]);
        viscous_force_y = global_visc_coeff*m_pd.laplacian(grids[mom_y]);
    }
    else if (m_viscosity_opt == "velocity"){
        viscous_force_x = global_visc_coeff*grids[rho]*m_pd.laplacian(grids[v_x]);
        viscous_force_y = global_visc_coeff*grids[rho]*m_pd.laplacian(grids[v_y]);
    }
    else{
        std::cerr << "Viscosity option <" << m_viscosity_opt << "> is not valid." << std::endl;
        assert(false);
    }
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


void IdealMHD2E::recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[i_press] = grids[n]*K_B*grids[i_temp].max(m_pd.temp_min);
    grids[e_press] = grids[n]*K_B*grids[e_temp].max(m_pd.temp_min);
    grids[i_thermal_energy] = (grids[i_press]/(m_pd.m_adiabatic_index - 1.0)).max(m_pd.thermal_energy_min);
    grids[e_thermal_energy] = (grids[e_press]/(m_pd.m_adiabatic_index - 1.0)).max(m_pd.thermal_energy_min);
}

void IdealMHD2E::recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[rho] = grids[n]*m_pd.m_ion_mass;
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[kinetic_energy] = 0.5*grids[rho]*(grids[v_x].square()+grids[v_y].square());
    grids[i_thermal_energy] = grids[i_thermal_energy].max(m_pd.thermal_energy_min);
    grids[e_thermal_energy] = grids[e_thermal_energy].max(m_pd.thermal_energy_min);
    grids[i_press] = (m_pd.m_adiabatic_index - 1.0)*grids[i_thermal_energy];
    grids[e_press] = (m_pd.m_adiabatic_index - 1.0)*grids[e_thermal_energy];
    grids[press] = grids[i_press] + grids[e_press];
    grids[i_temp] = (grids[i_press]/(K_B*grids[n])).max(m_pd.temp_min);
    grids[e_temp] = (grids[e_press]/(K_B*grids[n])).max(m_pd.temp_min);
    grids[b_x] = m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x];
    grids[b_y] = m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_mag] = (grids[b_x].square() + grids[b_y].square()).sqrt();
    grids[b_hat_x] = grids[b_x]/grids[b_mag];
    grids[b_hat_y] = grids[b_y]/grids[b_mag];
    grids[dPdx] = m_pd.derivative1D(grids[press],0);
    grids[dPdy] = m_pd.derivative1D(grids[press],1);
    catchNullFieldDirection(grids);
}

void IdealMHD2E::catchNullFieldDirection(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
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

    // m_grids[dt] = diagonals/(c_s + v_alfven + v_fast + v_slow + vel_mag);
    m_grids[dt] = diagonals/(vel_mag + c_s.max(v_alfven).max(v_fast).max(v_slow));
}

void IdealMHD2E::propagateChanges(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    m_pd.updateGhostZones();
    recomputeDerivedVarsFromEvolvedVars(grids);
    recomputeDT();
}

void IdealMHD2E::populateVariablesFromState(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeEvolvedVarsFromStateVars(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVarsFromEvolvedVars(grids);
    recomputeDT();
}

// sets viscosity in outermost interior cell to zero because it is large due to using a boundary cell
void IdealMHD2E::viscosity_mask(Grid& grid) const
{
    // treat left boundary
    for (int j=0; j<grid.cols(); j++){
        grid(m_pd.m_xl,j) = 0;
    }
    // treat right boundary
    for (int j=0; j<grid.cols(); j++){
        grid(m_pd.m_xu,j) = 0;
    }
    // treat bottom boundary
    for (int i=0; i<grid.rows(); i++){
        grid(i,m_pd.m_yl) = 0;
    }
    // treat top boundary
    for (int i=0; i<grid.rows(); i++){
        grid(i,m_pd.m_yu) = 0;
    }
}