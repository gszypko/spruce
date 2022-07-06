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

Grid IdealMHD::getDT(){
    return m_grids[dt];
}

void IdealMHD::propagateChanges(std::vector<Grid> &grids)
{
    catchUnderdensity(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVariables(grids);
}

void IdealMHD::populateVariablesFromState(std::vector<Grid> &grids){
    computeConstantTerms(grids);
    recomputeThermalEnergy(grids);
    recomputeDerivedVariables(grids);
    m_pd.updateGhostZones();
}

std::vector<Grid> IdealMHD::computeTimeDerivativesCharacteristicBoundary(const std::vector<Grid> &grids, bool x_bound_1, bool x_bound_2, bool y_bound_1, bool y_bound_2){
    std::vector<std::vector<Grid> > results(4, std::vector<Grid>(4,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim)));
    if(x_bound_1) results[0] = singleBoundaryTermsMOC(grids, 0, true);
    if(x_bound_2) results[1] = singleBoundaryTermsMOC(grids, 0, false);
    if(y_bound_1) results[2] = singleBoundaryTermsMOC(grids, 1, true);
    if(y_bound_2) results[3] = singleBoundaryTermsMOC(grids, 1, false);
    // for(Grid g : results[3]) std::cout << g << std::endl;
    // abort();
    Grid rho_sum = results[0][0] + results[1][0] + results[2][0] + results[3][0]; //all normal terms included
    rho_sum += results[0][1] + results[1][1] + results[2][1] + results[3][1]; //all parallel terms inside domain included

    // Only include corner parallel term in corners where only one characteristic boundary is active
    if(x_bound_1){
        if(!y_bound_1) rho_sum += results[0][2];
        if(!y_bound_2) rho_sum += results[0][3];
    }
    if(x_bound_2){
        if(!y_bound_1) rho_sum += results[1][2];
        if(!y_bound_2) rho_sum += results[1][3];
    }
    if(y_bound_1){
        if(!x_bound_1) rho_sum += results[2][2];
        if(!x_bound_2) rho_sum += results[2][3];
    }
    if(y_bound_2){
        if(!x_bound_1) rho_sum += results[3][2];
        if(!x_bound_2) rho_sum += results[3][3];
    }

    Grid zero = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    return {rho_sum, zero, zero, zero, zero, zero};

}

// Returns {normal term, parallel term inside domain, parallel term in lower corner, parallel term in upper corner} (in full-sized grid) for single MOC boundary
std::vector<Grid> IdealMHD::singleBoundaryTermsMOC(const std::vector<Grid> &grids, int boundary_index, bool boundary_lower){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    std::vector<int> x_idx, y_idx;
    if(boundary_index == 0){
        if(boundary_lower) for(int i=0; i<N_GHOST; i++) x_idx.push_back(i);
        else for(int i=N_GHOST; i>0; i--) x_idx.push_back(m_pd.m_xdim - i);
        for(int j=0; j<m_pd.m_ydim; j++) y_idx.push_back(j);
    } else {
        if(boundary_lower) for(int j=0; j<N_GHOST; j++) y_idx.push_back(j);
        else for(int j=N_GHOST; j>0; j--) y_idx.push_back(m_pd.m_ydim - j);
        for(int i=0; i<m_pd.m_xdim; i++) x_idx.push_back(i);
    }
    Grid boundaryMask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    for(int i : x_idx) for(int j : y_idx) boundaryMask(i,j) = 1.0;
    // for(int i : x_idx) std::cout << i << " ";
    // std::cout << std::endl;
    // for(int j : y_idx) std::cout << j << " ";
    // std::cout << std::endl;

    Grid parallel_term;
    if(boundary_index == 0) parallel_term = -m_pd.transportDerivative1D(grids[rho],grids[v_y],1,x_idx[0],y_idx[0]+N_GHOST,x_idx.back(),y_idx.back()-N_GHOST);
    else if(boundary_index == 1) parallel_term = -m_pd.transportDerivative1D(grids[rho],grids[v_x],0,x_idx[0]+N_GHOST,y_idx[0],x_idx.back()-N_GHOST,y_idx.back());
    else assert(false && "This function assumes two-dimensional grids");

    Grid parallel_term_corner_lower, parallel_term_corner_upper;
    if(boundary_index == 0){
        parallel_term_corner_lower = -m_pd.derivative1D(grids[mom_y],1,x_idx[0],y_idx[0]+1,x_idx.back(),y_idx[0]+1);
        parallel_term_corner_upper = -m_pd.derivative1D(grids[mom_y],1,x_idx[0],y_idx.back()-1,x_idx.back(),y_idx.back()-1);
    }
    else{
        parallel_term_corner_lower = -m_pd.derivative1D(grids[mom_x],0,x_idx[0]+1,y_idx[0],x_idx[0]+1,y_idx.back());
        parallel_term_corner_upper = -m_pd.derivative1D(grids[mom_x],0,x_idx.back()-1,y_idx[0],x_idx.back()-1,y_idx.back());
    }
    
    Grid normal_term;
    int v_perp = (boundary_index == 0) ? v_x : v_y;
    int v_para = (boundary_index == 1) ? v_x : v_y;
    bool positive_forward = !boundary_lower;

    Grid m_b_perp = (boundary_index == 0 ? (m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x]) : (m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y]));
    Grid m_b_para = (boundary_index == 1 ? (m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x]) : (m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y]));

    Grid s_perp = m_b_perp/m_b_perp.abs();
    Grid R_para = (m_b_para/m_b_para).abs();
    Grid c_s_sq = m_pd.m_adiabatic_index*grids[press]/grids[rho];
    Grid c_a_sq = m_grids[b_mag].square()/(4.0*PI*m_grids[rho]);
    Grid c_perp_sq = m_b_perp.square()/(4.0*PI*m_grids[rho]);
    Grid c_plus_sq = (c_a_sq + c_s_sq)/2.0 + ((c_a_sq + c_s_sq).square()/4.0 - c_perp_sq*c_s_sq).sqrt();
    Grid c_minus_sq = (c_a_sq + c_s_sq)/2.0 - ((c_a_sq + c_s_sq).square()/4.0 - c_perp_sq*c_s_sq).sqrt();
    Grid alpha_plus_sq = ((c_s_sq - c_minus_sq)/(c_plus_sq - c_minus_sq)).max(0.0);
    Grid alpha_minus_sq = ((c_plus_sq - c_s_sq)/(c_plus_sq - c_minus_sq)).max(0.0);
    Grid c_s = c_s_sq.sqrt();
    Grid c_plus = c_plus_sq.sqrt(), c_minus = c_minus_sq.sqrt();
    Grid alpha_plus = alpha_plus_sq.sqrt(), alpha_minus = alpha_minus_sq.sqrt();
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(m_b_perp(i,j) == 0.0) s_perp(i,j) = 0.0;
            if(m_b_para(i,j) == 0.0) R_para(i,j) = 0.0;
        }
    }
    // std::cout << "c_plus\n";
    // std::cout << alpha_plus_sq << std::endl;
    // std::cout << alpha_minus_sq << std::endl;
    // std::cout << alpha_plus << std::endl;
    // std::cout << alpha_minus << std::endl;
    Grid d_2 = grids[v_perp]*(
        (grids[thermal_energy] + grids[press])*m_pd.derivative1DBackward(grids[rho],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - grids[rho]*m_pd.derivative1DBackward(grids[thermal_energy],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
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
    // if(m_pd.m_iter == 10){
    //     std::cout << "before\n";
    //     std::cout << d_2 << std::endl;
    //     std::cout << d_5 << std::endl;
    //     std::cout << d_6 << std::endl;
    //     std::cout << d_7 << std::endl;
    //     std::cout << d_8 << std::endl;
    // }
    //Compute remaining z-derivative terms
    int inflow_dir = boundary_lower ? 1 : -1;
    for(int i : x_idx){
        for(int j : y_idx){
            //Sift through to negate any d terms for which the eigenvalue is inward-pointing (REMEMBER TO CHANGE THIS BEHAVIOR IN THE FUTURE FOR OTHER EQUATIONS)
            if(inflow_dir*grids[v_perp](i,j) >= 0.0) d_2(i,j) = 0.0;
            if(inflow_dir*(grids[v_perp](i,j) + c_plus(i,j)) >= 0.0) d_5(i,j) = 0.0;
            if(inflow_dir*(grids[v_perp](i,j) - c_plus(i,j)) >= 0.0) d_6(i,j) = 0.0;
            if(inflow_dir*(grids[v_perp](i,j) + c_minus(i,j)) >= 0.0) d_7(i,j) = 0.0;
            if(inflow_dir*(grids[v_perp](i,j) - c_minus(i,j)) >= 0.0) d_8(i,j) = 0.0;
        }
    }
    normal_term = -( (m_pd.m_adiabatic_index/grids[rho])*d_2 + 0.5*grids[rho]*alpha_plus*(d_5+d_6) + 0.5*grids[rho]*alpha_minus*(d_7+d_8) )/c_s_sq;

    return {normal_term,parallel_term,parallel_term_corner_lower,parallel_term_corner_upper};
}