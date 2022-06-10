#include "idealmhd.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"

IdealMHD::IdealMHD(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void IdealMHD::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[rho] += step*m_pd.m_ghost_zone_mask*time_derivatives[0];
    grids[mom_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[1];
    grids[mom_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[2];
    grids[bi_x] += step*m_pd.m_ghost_zone_mask*time_derivatives[3];
    grids[bi_y] += step*m_pd.m_ghost_zone_mask*time_derivatives[4];
    grids[thermal_energy] += step*m_pd.m_ghost_zone_mask*time_derivatives[5];
    propagateChanges(grids);
}

std::vector<Grid> IdealMHD::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    //Advance time by min_dt
    const Grid &m_mom_x = grids[mom_x], &m_mom_y = grids[mom_y],
        &m_v_x = grids[v_x], &m_v_y = grids[v_y],
        &m_be_x = m_pd.m_internal_grids[PlasmaDomain::be_x], &m_be_y = m_pd.m_internal_grids[PlasmaDomain::be_y],
        &m_bi_x = grids[bi_x], &m_bi_y = grids[bi_y],
        &m_rho = grids[rho], &m_thermal_energy = grids[thermal_energy], &m_press = grids[press];

    std::vector<Grid> m_v = {m_v_x,m_v_y};

    Grid viscous_force_x = visc_coeff*m_pd.laplacian(m_mom_x);
    Grid viscous_force_y = visc_coeff*m_pd.laplacian(m_mom_y);

    Grid d_rho_dt = -m_pd.transportDivergence2D(m_rho,m_v);

    Grid curl_db = m_pd.curl2D(m_bi_x,m_bi_y)/(4.0*PI);
    std::vector<Grid> mag_mom_terms_firstorder = Grid::CrossProductZ2D(curl_db,{m_be_x,m_be_y});
    std::vector<Grid> mag_mom_terms_secorder = Grid::CrossProductZ2D(curl_db,{m_bi_x,m_bi_y}); //second order terms
    Grid d_mom_x_dt = -m_pd.transportDivergence2D(m_mom_x, m_v)
                    - m_pd.derivative1D(m_press, 0)
                    + m_pd.m_ghost_zone_mask * (m_rho*grids[grav_x] + viscous_force_x
                    + mag_mom_terms_firstorder[0] + mag_mom_terms_secorder[0]);
    Grid d_mom_y_dt = -m_pd.transportDivergence2D(m_mom_y, m_v)
                    - m_pd.derivative1D(m_press, 1)
                    + m_pd.m_ghost_zone_mask * (m_rho*grids[grav_y] + viscous_force_y
                    + mag_mom_terms_firstorder[1] + mag_mom_terms_secorder[1]);

    std::vector<Grid> induction_rhs_b0 = m_pd.curlZ(Grid::CrossProduct2D({m_v_x,m_v_y},{m_be_x,m_be_y}));
    std::vector<Grid> induction_rhs_db = m_pd.curlZ(Grid::CrossProduct2D({m_v_x,m_v_y},{m_bi_x,m_bi_y}));
    Grid d_bi_x_dt = induction_rhs_b0[0] + induction_rhs_db[0];
    Grid d_bi_y_dt = induction_rhs_b0[1] + induction_rhs_db[1];

    Grid d_thermal_energy_dt = - m_pd.transportDivergence2D(m_thermal_energy,m_v)
                                - m_press*m_pd.divergence2D(m_v);

    return {d_rho_dt, d_mom_x_dt, d_mom_y_dt, d_bi_x_dt, d_bi_y_dt, d_thermal_energy_dt};
}

void IdealMHD::computeConstantTerms(std::vector<Grid> &grids){
    grids[div_be] = m_pd.divergence2D(m_pd.m_internal_grids[PlasmaDomain::be_x],m_pd.m_internal_grids[PlasmaDomain::be_y]);
}

void IdealMHD::recomputeDerivedVariables(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[press] = 2.0*K_B*grids[rho]*grids[temp]/m_pd.m_ion_mass;
    grids[thermal_energy] = grids[press]/(m_pd.m_adiabatic_index - 1.0);
    grids[kinetic_energy] = 0.5*(grids[mom_x].square() + grids[mom_y].square())/grids[rho];
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[n] = grids[rho]/m_pd.m_ion_mass;
    grids[div_bi] = m_pd.divergence2D(grids[bi_x],grids[bi_y]);
    Grid m_b_x = m_pd.m_internal_grids[PlasmaDomain::be_x] + grids[bi_x], m_b_y = m_pd.m_internal_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_magnitude] = (m_b_x.square() + m_b_y.square()).sqrt();
    //Need to ensure that b_hat is zero when b is zero
    Grid &b_h_x = grids[b_hat_x], &b_h_y = grids[b_hat_y], &b_mag = grids[b_magnitude];
    b_h_x = m_b_x/b_mag;
    b_h_y = m_b_y/b_mag;
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(b_mag(i,j) == 0.0){
                b_h_x(i,j) = 0.0;
                b_h_y(i,j) = 0.0;
            }
        }
    }
    recomputeDT();
}

void IdealMHD::recomputeTemperature(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    Grid &m_press = grids[press], &m_thermal_energy = grids[thermal_energy];
    m_thermal_energy = m_thermal_energy.max(m_pd.thermal_energy_min);
    m_press = (m_pd.m_adiabatic_index - 1.0)*m_thermal_energy;
    grids[temp] = (m_pd.m_ion_mass*m_press/(2.0*K_B*grids[rho])).max(m_pd.temp_min);
} //need to rethink this in the general case

void IdealMHD::catchUnderdensity(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(grids[rho](i,j) < m_pd.rho_min){
                if(i >= m_pd.m_xl && i <= m_pd.m_xu && j >= m_pd.m_yl && j <= m_pd.m_yu){
                    grids[rho](i,j) = m_pd.rho_min;
                }
            }
        }
    }
}

void IdealMHD::recomputeDT(){
    Grid c_s = (m_pd.m_adiabatic_index*m_grids[press]/m_grids[rho]).sqrt();
    Grid c_s_sq = c_s.square();

    Grid v_alfven = m_grids[b_magnitude]/(4.0*PI*m_grids[rho]).sqrt();
    m_grids[v_a] = v_alfven;
    Grid v_alfven_sq = v_alfven.square();
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);

    //Magnetoacoustic modes
    Grid delta = (one - 4.0*c_s_sq*v_alfven_sq/(c_s_sq+v_alfven_sq).square()).sqrt();
    Grid v_fast = (0.5*(c_s_sq + v_alfven_sq)*(one + delta)).sqrt();
    Grid v_slow = (0.5*(c_s_sq + v_alfven_sq)*(one - delta)).sqrt();

    //Bulk velocity transit time
    Grid diagonals = (m_pd.m_internal_grids[PlasmaDomain::d_x].square() + m_pd.m_internal_grids[PlasmaDomain::d_y].square()).sqrt();
    Grid vel_mag = (m_grids[v_x].square() + m_grids[v_y].square()).sqrt();

    m_grids[dt] = diagonals/(c_s + v_alfven + v_fast + v_slow + vel_mag);
}

Grid IdealMHD::getDT(){
    return m_grids[dt];
}

void IdealMHD::propagateChanges(std::vector<Grid> &grids)
{
    catchUnderdensity(grids);
    m_pd.updateGhostZones();
    recomputeTemperature(grids);
    recomputeDerivedVariables(grids);
}

void IdealMHD::populateVariablesFromState(std::vector<Grid> &grids){
    computeConstantTerms(grids);
    recomputeDerivedVariables(grids);
    m_pd.updateGhostZones();
}
