#include "idealmhd.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>

IdealMHD::IdealMHD(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void IdealMHD::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    std::vector<Grid> char_evolution = computeTimeDerivativesCharacteristicBoundary(grids,
                                                m_pd.x_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC,
                                                m_pd.x_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC,
                                                m_pd.y_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC,
                                                m_pd.y_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC);
    grids[rho] += step*(m_pd.m_ghost_zone_mask*time_derivatives[0] + char_evolution[0]);
    grids[mom_x] += step*(m_pd.m_ghost_zone_mask*time_derivatives[1] + char_evolution[1]);
    grids[mom_y] += step*(m_pd.m_ghost_zone_mask*time_derivatives[2] + char_evolution[2]);
    grids[thermal_energy] += step*(m_pd.m_ghost_zone_mask*time_derivatives[3] + char_evolution[3]);
    grids[bi_x] += step*(m_pd.m_ghost_zone_mask*time_derivatives[4] + char_evolution[4]);
    grids[bi_y] += step*(m_pd.m_ghost_zone_mask*time_derivatives[5] + char_evolution[5]);
    propagateChanges(grids);
}

std::vector<Grid> IdealMHD::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff){
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
    std::vector<Grid> internal_mag_force = Grid::CrossProductZ2D(curl_db,{grids[bi_x],grids[bi_y]});
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
    Grid d_thermal_energy_dt =  - m_pd.transportDivergence2D(grids[thermal_energy],v)
                                - grids[press]*m_pd.divergence2D(v);
    // induction equations
    std::vector<Grid> induction_rhs_external = m_pd.curlZ(Grid::CrossProduct2D(v,{be_x,be_y}));
    std::vector<Grid> induction_rhs_internal = m_pd.curlZ(Grid::CrossProduct2D(v,{grids[bi_x],grids[bi_y]}));
    Grid d_bi_x_dt = induction_rhs_external[0] + induction_rhs_internal[0];
    Grid d_bi_y_dt = induction_rhs_external[1] + induction_rhs_internal[1];
    // return time derivatives
    return {d_rho_dt, d_mom_x_dt, d_mom_y_dt, d_thermal_energy_dt, d_bi_x_dt, d_bi_y_dt};
}

void IdealMHD::computeConstantTerms(std::vector<Grid> &grids){
}

void IdealMHD::recomputeDerivedVariables(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    recomputeTemperature(grids);
    recomputeKineticEnergy(grids);
    recomputeMagneticFields(grids);
    recomputeDT();
}

void IdealMHD::recomputeTemperature(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    // enforce minimum on thermal energies
    grids[thermal_energy] = grids[thermal_energy].max(m_pd.thermal_energy_min);
    // recompute pressures
    grids[press] = (m_pd.m_adiabatic_index - 1.0)*grids[thermal_energy];
    // recompute temperatures and enforce minimum
    grids[temp] = (grids[press]/(2*K_B*grids[n])).max(m_pd.temp_min);
}

void IdealMHD::recomputeThermalEnergy(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeNumberDensity(grids);
    grids[press] = 2*K_B*grids[n]*grids[temp];
    grids[thermal_energy] = grids[press]/(m_pd.m_adiabatic_index - 1.0);
}

void IdealMHD::recomputeKineticEnergy(std::vector<Grid> &grids){
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[kinetic_energy] = 0.5*grids[rho]*(grids[v_x].square()+grids[v_y].square());
}

void IdealMHD::recomputeNumberDensity(std::vector<Grid> &grids){
    grids[n] = grids[rho]/m_pd.m_ion_mass;
}

void IdealMHD::recomputeMagneticFields(std::vector<Grid> &grids){
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

void IdealMHD::recomputeDT(){
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

Grid IdealMHD::getDT(){
    return m_grids[dt];
}

void IdealMHD::enforceMinimums(std::vector<Grid>& grids)
{
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[rho] = grids[n]*m_pd.m_ion_mass;
    grids[thermal_energy] = grids[thermal_energy].max(m_pd.thermal_energy_min);
}

void IdealMHD::propagateChanges(std::vector<Grid> &grids)
{
    enforceMinimums(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVariables(grids);
}

void IdealMHD::populateVariablesFromState(std::vector<Grid> &grids){
    computeConstantTerms(grids);
    recomputeThermalEnergy(grids);
    recomputeDerivedVariables(grids);
    m_pd.updateGhostZones();
}

// Returns time derivative from characteristic boundary cond for the quantities
// {rho, mom_x, mom_y, thermal_energy, b_x, b_y} in that order
std::vector<Grid> IdealMHD::computeTimeDerivativesCharacteristicBoundary(const std::vector<Grid> &grids, bool x_bound_1, bool x_bound_2, bool y_bound_1, bool y_bound_2){
    std::vector<std::vector<Grid> > x_lower, x_upper, y_lower, y_upper;
    if(x_bound_1) x_lower = singleBoundaryTermsMOC(grids, 0, true);
    if(x_bound_2) x_upper = singleBoundaryTermsMOC(grids, 0, false);
    if(y_bound_1) y_lower = singleBoundaryTermsMOC(grids, 1, true);
    if(y_bound_2) y_upper = singleBoundaryTermsMOC(grids, 1, false);
    int num_vars = x_lower.size();
    std::vector<Grid> results(num_vars,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    for(int i=0; i<num_vars; i++){
        results[i] = x_lower[i][0] + x_upper[i][0] + y_lower[i][0] + y_upper[i][0]; //all normal terms included
        results[i] += x_lower[i][1] + x_upper[i][1] + y_lower[i][1] + y_upper[i][1]; //all parallel terms inside domain included
    }

    // Only include corner parallel term in corners where only one characteristic boundary is active
    std::vector<int> vars = {0,1,2,3,4,5};
    if(x_bound_1){
        if(!y_bound_1) for(int i : vars) results[i] += x_lower[i][2];
        if(!y_bound_2) for(int i : vars) results[i] += x_lower[i][3];
    }
    if(x_bound_2){
        if(!y_bound_1) for(int i : vars) results[i] += x_upper[i][2];
        if(!y_bound_2) for(int i : vars) results[i] += x_upper[i][3];
    }
    if(y_bound_1){
        if(!x_bound_1) for(int i : vars) results[i] += y_lower[i][2];
        if(!x_bound_2) for(int i : vars) results[i] += y_lower[i][3];
    }
    if(y_bound_2){
        if(!x_bound_1) for(int i : vars) results[i] += y_upper[i][2];
        if(!x_bound_2) for(int i : vars) results[i] += y_upper[i][3];
    }

    Grid mom_x_term = grids[rho]*results[1] + grids[v_x]*results[0];
    Grid mom_y_term = grids[rho]*results[2] + grids[v_y]*results[0];
    results[1] = mom_x_term;
    results[2] = mom_y_term;

    return results;
}

// Returns {normal term, parallel term inside domain, parallel term in lower corner, parallel term in upper corner} (in full-sized grid) for single MOC boundary,
// corresponding to each of {rho, v_x, v_y, b_x, b_y, thermal_energy}
std::vector<std::vector<Grid>> IdealMHD::singleBoundaryTermsMOC(const std::vector<Grid> &grids, int boundary_index, bool boundary_lower){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    std::vector<int> x_idx, y_idx; //Full extent of boundary's ghost zone region
    std::vector<int> x_idx_interior, y_idx_interior; //Full extent of boundary's ghost zone not overlapping with other ghost zones
    std::vector<int> x_idx_padded, y_idx_padded; //Full extent of boundary's ghost zone region, reduced for central differencing safety
    std::vector<int> x_idx_lower, y_idx_lower; //Lower corner of boundary's ghost zone, central difference safe
    std::vector<int> x_idx_upper, y_idx_upper; //Upper corner of boundary's ghost zone, central difference safe
    if(boundary_index == 0){
        if(boundary_lower) x_idx = {0, N_GHOST-1};
        else x_idx = {(int)(m_pd.m_xdim) - N_GHOST, (int)(m_pd.m_xdim) - 1};
        y_idx = {0, (int)(m_pd.m_ydim) - 1};
    } else {
        if(boundary_lower) y_idx = {0, N_GHOST-1};
        else y_idx = {(int)(m_pd.m_ydim) - N_GHOST, (int)(m_pd.m_ydim) - 1};
        x_idx = {0, (int)(m_pd.m_xdim) - 1};
    }
    x_idx_interior = x_idx;
    x_idx_padded = x_idx;
    x_idx_lower = {0,-1};
    x_idx_upper = {0,-1};
    y_idx_interior = y_idx;
    y_idx_padded = y_idx;
    y_idx_lower = {0,-1};
    y_idx_upper = {0,-1};
    if(boundary_index == 0 && m_pd.y_bound_1 != PlasmaDomain::BoundaryCondition::Periodic){
        y_idx_interior = {y_idx[0]+N_GHOST,y_idx.back()-N_GHOST};
        y_idx_padded = {y_idx[0]+1,y_idx.back()-1};
        y_idx_lower = {y_idx[0],y_idx[0]+N_GHOST-1};
        y_idx_upper = {y_idx.back()-N_GHOST+1,y_idx.back()};
    } else if(boundary_index == 1 && m_pd.x_bound_1 != PlasmaDomain::BoundaryCondition::Periodic){
        x_idx_interior = {x_idx[0]+N_GHOST,x_idx.back()-N_GHOST};
        x_idx_padded = {x_idx[0]+1,x_idx.back()-1};
        x_idx_lower = {x_idx[0],x_idx[0]+N_GHOST-1};
        x_idx_upper = {x_idx.back()-N_GHOST+1,x_idx.back()};
    }

    // Mask grid where equal to 1 inside of boundary ghost zone, 0 elsewhere
    Grid boundary_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid interior_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid lower_corner_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid upper_corner_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    for(int i=x_idx[0]; i<=x_idx[1]; i++) for(int j=y_idx[0]; j<=y_idx[1]; j++) boundary_mask(i,j) = 1.0;
    for(int i=x_idx_interior[0]; i<=x_idx_interior[1]; i++) for(int j=y_idx_interior[0]; j<=y_idx_interior[1]; j++) interior_mask(i,j) = 1.0;
    for(int i=x_idx_lower[0]; i<=x_idx_lower[1]; i++) for(int j=y_idx_lower[0]; j<=y_idx_lower[1]; j++) lower_corner_mask(i,j) = 1.0;
    for(int i=x_idx_upper[0]; i<=x_idx_upper[1]; i++) for(int j=y_idx_upper[0]; j<=y_idx_upper[1]; j++) upper_corner_mask(i,j) = 1.0;

    int v_perp = (boundary_index == 0) ? v_x : v_y;
    int v_para = (boundary_index == 1) ? v_x : v_y;
    int grav_perp = (boundary_index == 0) ? grav_x : grav_y;
    int grav_para = (boundary_index == 1) ? grav_x : grav_y;
    bool positive_forward = !boundary_lower;

    Grid m_b_perp = (boundary_index == 0 ? (m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x]) : (m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y]));
    Grid m_b_para = (boundary_index == 1 ? (m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x]) : (m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y]));

    Grid s_perp = m_b_perp/(m_b_perp.abs());
    Grid R_para = m_b_para/(m_b_para.abs());
    Grid c_s_sq = m_pd.m_adiabatic_index*grids[press]/grids[rho];
    Grid c_a_sq = m_grids[b_mag].square()/(4.0*PI*m_grids[rho]);
    Grid c_perp_sq = m_b_perp.square()/(4.0*PI*m_grids[rho]);
    Grid c_plus_sq = (c_a_sq + c_s_sq)/2.0 + ((c_a_sq + c_s_sq).square()/4.0 - c_perp_sq*c_s_sq).sqrt();
    Grid c_minus_sq = (c_a_sq + c_s_sq)/2.0 - ((c_a_sq + c_s_sq).square()/4.0 - c_perp_sq*c_s_sq).sqrt();
    Grid alpha_plus_sq = ((c_s_sq - c_minus_sq)/(c_plus_sq - c_minus_sq)).max(0.0);
    Grid alpha_minus_sq = ((c_plus_sq - c_s_sq)/(c_plus_sq - c_minus_sq)).max(0.0);
    Grid c_s = c_s_sq.sqrt();
    Grid c_perp = c_perp_sq.sqrt();
    Grid c_plus = c_plus_sq.sqrt(), c_minus = c_minus_sq.sqrt();
    Grid alpha_plus = alpha_plus_sq.sqrt(), alpha_minus = alpha_minus_sq.sqrt();
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(m_b_perp(i,j) == 0.0) s_perp(i,j) = 0.0;
            if(m_b_para(i,j) == 0.0) R_para(i,j) = 1.0;
        }
    }

    //Compute ordinary perpendicular terms
    Grid d_1 = grids[v_perp]*(
        m_pd.derivative1DBackward(m_b_perp,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );
    Grid d_2 = grids[v_perp]*(
        (grids[thermal_energy] + grids[press])*m_pd.derivative1DBackward(grids[rho],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - grids[rho]*m_pd.derivative1DBackward(grids[thermal_energy],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );
    Grid d_3 = (grids[v_perp] + c_perp)*(
        Grid::Zero(m_pd.m_xdim,m_pd.m_ydim)
    );
    Grid d_4 = (grids[v_perp] - c_perp)*(
        Grid::Zero(m_pd.m_xdim,m_pd.m_ydim)
    );
    Grid d_5 = (grids[v_perp] + c_plus)*(
        alpha_plus/grids[rho]*m_pd.derivative1DBackward(grids[press],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - s_perp*R_para*c_minus*alpha_minus*m_pd.derivative1DBackward(grids[v_para],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        + c_plus*alpha_plus*m_pd.derivative1DBackward(grids[v_perp],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        + R_para*c_s*alpha_minus/(4.0*PI*(grids[rho])).sqrt()*m_pd.derivative1DBackward(m_b_para,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );
    Grid d_6 = (grids[v_perp] - c_plus)*(
        alpha_plus/grids[rho]*m_pd.derivative1DBackward(grids[press],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        + s_perp*R_para*c_minus*alpha_minus*m_pd.derivative1DBackward(grids[v_para],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - c_plus*alpha_plus*m_pd.derivative1DBackward(grids[v_perp],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        + R_para*c_s*alpha_minus/(4.0*PI*(grids[rho])).sqrt()*m_pd.derivative1DBackward(m_b_para,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );
    Grid d_7 = (grids[v_perp] + c_minus)*(
        alpha_minus/grids[rho]*m_pd.derivative1DBackward(grids[press],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        + s_perp*R_para*c_plus*alpha_plus*m_pd.derivative1DBackward(grids[v_para],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        + c_minus*alpha_minus*m_pd.derivative1DBackward(grids[v_perp],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - R_para*c_s*alpha_plus/(4.0*PI*(grids[rho])).sqrt()*m_pd.derivative1DBackward(m_b_para,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );
    Grid d_8 = (grids[v_perp] - c_minus)*(
        alpha_minus/grids[rho]*m_pd.derivative1DBackward(grids[press],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - s_perp*R_para*c_plus*alpha_plus*m_pd.derivative1DBackward(grids[v_para],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - c_minus*alpha_minus*m_pd.derivative1DBackward(grids[v_perp],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - R_para*c_s*alpha_plus/(4.0*PI*(grids[rho])).sqrt()*m_pd.derivative1DBackward(m_b_para,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );

    int para_index = boundary_index == 0 ? 1 : 0;
    //Compute alternate perpendicular terms (for negating inflowing characteristics)
    Grid d_3_alt = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid d_5_alt = c_plus*alpha_plus*(
        grids[grav_perp]
        + 1.0/(4.0*PI*grids[rho])* m_b_para * m_pd.derivative1D(m_b_perp,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1]))
        + c_minus*alpha_minus*s_perp*(
            R_para * grids[grav_para]
            + R_para / grids[rho] * m_pd.derivative1D(grids[press] + m_b_perp.square()/(8.0*PI),para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1]));
    Grid d_7_alt = c_minus*alpha_minus*(
        grids[grav_perp]
        + 1.0/(4.0*PI*grids[rho])* m_b_para * m_pd.derivative1D(m_b_perp,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1]))
        - c_plus*alpha_plus*s_perp*(
            R_para * grids[grav_para]
            + R_para / grids[rho] * m_pd.derivative1D(grids[press] + m_b_perp.square()/(8.0*PI),para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1]));

    //Check for inflowing characteristic speeds and swap out alternate values as necessary
    int inflow_dir = boundary_lower ? 1 : -1;
    for(int i = x_idx[0]; i <= x_idx[1]; i++){
        for(int j = y_idx[0]; j <= y_idx[1]; j++){
            if(inflow_dir*grids[v_perp](i,j) > 0.0) { d_1(i,j) = 0.0; d_2(i,j) = 0.0; }
            if(inflow_dir*(grids[v_perp](i,j) + c_perp(i,j)) > 0.0) {d_3(i,j) = d_3_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) - c_perp(i,j)) > 0.0) {d_4(i,j) = -d_3_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) + c_plus(i,j)) > 0.0) {d_5(i,j) = d_5_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) - c_plus(i,j)) > 0.0) {d_6(i,j) = -d_5_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) + c_minus(i,j)) > 0.0) {d_7(i,j) = d_7_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) - c_minus(i,j)) > 0.0) {d_8(i,j) = -d_7_alt(i,j);}
        }
    }
    //{rho, v_x, v_y, b_x, b_y, thermal_energy}
    std::vector<Grid> rho_terms, v_para_terms, v_perp_terms, b_perp_terms, b_para_terms, thermal_energy_terms;

    // Compute all normal terms

    double visc_coeff = 0.0;
    rho_terms.push_back(-((m_pd.m_adiabatic_index/grids[rho])*d_2 + 0.5*grids[rho]*alpha_plus*(d_5+d_6) + 0.5*grids[rho]*alpha_minus*(d_7+d_8))/c_s_sq);
    v_para_terms.push_back(-0.5*s_perp*(c_minus*alpha_minus/c_s_sq*R_para*(-d_5 + d_6) + c_plus*alpha_plus/c_s_sq*R_para*(d_7 - d_8))
                            + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_para],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
                            + visc_coeff/grids[rho]*grids[v_para]*m_pd.secondDerivative1DBackward(grids[rho],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back()));
    v_perp_terms.push_back(-0.5/c_s_sq*(c_plus*alpha_plus*(d_5 - d_6) + c_minus*alpha_minus*(d_7 - d_8))
                            + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_perp],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
                            + visc_coeff/grids[rho]*grids[v_perp]*m_pd.secondDerivative1DBackward(grids[rho],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back()));
    thermal_energy_terms.push_back(-(0.5*(grids[thermal_energy]+grids[press])*alpha_plus*(d_5+d_6) + 0.5*(grids[thermal_energy]+grids[press])*alpha_minus*(d_7+d_8))/c_s_sq);
    b_para_terms.push_back(-(PI*grids[rho]).sqrt()*( alpha_minus/c_s*R_para*(d_5 + d_6) - alpha_plus/c_s*R_para*(d_7 + d_8) ));
    b_perp_terms.push_back(-d_1);


    // Compute all interior parallel terms
    rho_terms.push_back(-m_pd.derivative1D(grids[rho]*grids[v_para], para_index,
                        x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    v_para_terms.push_back(-1.0/grids[rho]*m_pd.derivative1D(grids[press] + (m_b_para.square() + m_b_perp.square())/(8.0*PI),para_index,
                                                            x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                        -grids[v_para]*m_pd.derivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                        +m_b_para/(4.0*PI*grids[rho])*m_pd.derivative1D(m_b_para,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                        +interior_mask*grids[grav_para]
                        + visc_coeff*m_pd.secondDerivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                        + visc_coeff/grids[rho]*grids[v_para]*m_pd.secondDerivative1D(grids[rho],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    v_perp_terms.push_back(-grids[v_para]*m_pd.derivative1D(grids[v_perp],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                        +m_b_para/(4.0*PI*grids[rho])*m_pd.derivative1D(m_b_perp,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                        +interior_mask*grids[grav_perp]
                        + visc_coeff*m_pd.secondDerivative1D(grids[v_perp],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                        + visc_coeff/grids[rho]*grids[v_perp]*m_pd.secondDerivative1D(grids[rho],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    thermal_energy_terms.push_back(-m_pd.derivative1D(grids[thermal_energy]*grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                -grids[press]*m_pd.derivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    b_para_terms.push_back(-grids[v_para]*m_pd.derivative1D(m_b_para,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                            -m_b_para*m_pd.derivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                            +m_b_para*m_pd.derivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    b_perp_terms.push_back(-grids[v_para]*m_pd.derivative1D(m_b_perp,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                            -m_b_perp*m_pd.derivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                            +m_b_para*m_pd.derivative1D(grids[v_perp],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));

    // Compute all lower corner parallel terms

    rho_terms.push_back(-m_pd.derivative1DBackward(grids[boundary_index == 0 ? mom_y : mom_x], false, para_index,
                        x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1]));
    v_para_terms.push_back(-1.0/grids[rho]*m_pd.derivative1DBackward(grids[press] + (m_b_para.square() + m_b_perp.square())/(8.0*PI),false, para_index,
                                                            x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                        -grids[v_para]*m_pd.derivative1DBackward(grids[v_para],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                        +m_b_para/(4.0*PI*grids[rho])*m_pd.derivative1DBackward(m_b_para,false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                        +lower_corner_mask*grids[grav_para]
                        + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_para],false,para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                        + visc_coeff/grids[rho]*grids[v_para]*m_pd.secondDerivative1DBackward(grids[rho],false,para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1]));
    v_perp_terms.push_back(-grids[v_para]*m_pd.derivative1DBackward(grids[v_perp],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                        +m_b_para/(4.0*PI*grids[rho])*m_pd.derivative1DBackward(m_b_perp,false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                        +lower_corner_mask*grids[grav_perp]
                        + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_perp],false,para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                        + visc_coeff/grids[rho]*grids[v_perp]*m_pd.secondDerivative1DBackward(grids[rho],false,para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1]));
    thermal_energy_terms.push_back(-m_pd.derivative1DBackward(grids[thermal_energy]*grids[v_para],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                                -grids[press]*m_pd.derivative1DBackward(grids[v_para],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1]));
    b_para_terms.push_back(-grids[v_para]*m_pd.derivative1DBackward(m_b_para,false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                            -m_b_para*m_pd.derivative1DBackward(grids[v_para],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                            +m_b_para*m_pd.derivative1DBackward(grids[v_para],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1]));
    b_perp_terms.push_back(-grids[v_para]*m_pd.derivative1DBackward(m_b_perp,false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                            -m_b_perp*m_pd.derivative1DBackward(grids[v_para],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1])
                            +m_b_para*m_pd.derivative1DBackward(grids[v_perp],false, para_index,x_idx_lower[0],y_idx_lower[0],x_idx_lower[1],y_idx_lower[1]));

    // Compute all upper corner parallel terms

    rho_terms.push_back(-m_pd.derivative1DBackward(grids[boundary_index == 0 ? mom_y : mom_x], true, para_index,
                        x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1]));
    v_para_terms.push_back(-1.0/grids[rho]*m_pd.derivative1DBackward(grids[press] + (m_b_para.square() + m_b_perp.square())/(8.0*PI),true, para_index,
                                                            x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                        -grids[v_para]*m_pd.derivative1DBackward(grids[v_para],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                        +m_b_para/(4.0*PI*grids[rho])*m_pd.derivative1DBackward(m_b_para,true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                        +upper_corner_mask*grids[grav_para]
                        + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_para],true,para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                        + visc_coeff/grids[rho]*grids[v_para]*m_pd.secondDerivative1DBackward(grids[rho],true,para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1]));
    v_perp_terms.push_back(-grids[v_para]*m_pd.derivative1DBackward(grids[v_perp],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                        +m_b_para/(4.0*PI*grids[rho])*m_pd.derivative1DBackward(m_b_perp,true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                        +upper_corner_mask*grids[grav_perp]
                        + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_perp],true,para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                        + visc_coeff/grids[rho]*grids[v_perp]*m_pd.secondDerivative1DBackward(grids[rho],true,para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1]));
    thermal_energy_terms.push_back(-m_pd.derivative1DBackward(grids[thermal_energy]*grids[v_para],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                                -grids[press]*m_pd.derivative1DBackward(grids[v_para],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1]));
    b_para_terms.push_back(-grids[v_para]*m_pd.derivative1DBackward(m_b_para,true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                            -m_b_para*m_pd.derivative1DBackward(grids[v_para],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                            +m_b_para*m_pd.derivative1DBackward(grids[v_para],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1]));
    b_perp_terms.push_back(-grids[v_para]*m_pd.derivative1DBackward(m_b_perp,true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                            -m_b_perp*m_pd.derivative1DBackward(grids[v_para],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1])
                            +m_b_para*m_pd.derivative1DBackward(grids[v_perp],true, para_index,x_idx_upper[0],y_idx_upper[0],x_idx_upper[1],y_idx_upper[1]));

    std::vector<Grid> zero(4,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    return {rho_terms,
            v_x == v_para ? v_para_terms : v_perp_terms, v_y == v_para ? v_para_terms : v_perp_terms,
            thermal_energy_terms,
            boundary_index == 0 ? b_perp_terms : b_para_terms, boundary_index == 1 ? b_perp_terms : b_para_terms};
}

