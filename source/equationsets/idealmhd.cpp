#include "idealmhd.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>
#include <cmath>
#include <cassert>

IdealMHD::IdealMHD(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void IdealMHD::parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) 
{
    for (int i=0; i<lhs.size(); i++){
        if (lhs[i] == "global_viscosity") m_global_viscosity = stod(rhs[i]);
        else if (lhs[i] == "viscosity_opt") m_viscosity_opt = rhs[i];
        else if (lhs[i] == "moc_b_limiting") moc_b_limiting = (rhs[i] == "true");
        else if (lhs[i] == "moc_mom_limiting") moc_mom_limiting = (rhs[i] == "true");
        else if (lhs[i] == "moc_b_lim_threshold"){
            moc_b_lim_threshold = stod(rhs[i]);
            assert(moc_b_lim_threshold <= 1.0 && "MoC B field threshold should be <=1.0");
        }
        else if (lhs[i] == "moc_mom_lim_threshold"){
            moc_mom_lim_threshold = stod(rhs[i]);
            assert(moc_mom_lim_threshold <= 1.0 && "MoC momentum threshold should be <=1.0");
        }
        else{
            std::cerr << lhs[i] << " is not recognized for this equation set." << std::endl;
            assert(false);
        }
    }
}

void IdealMHD::setupEquationSetDerived()
{
}

std::vector<Grid> IdealMHD::computeTimeDerivativesDerived(const std::vector<Grid> &grids) const {
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // PlasmaDomain grid references for more concise notation
    const Grid& d_x = m_pd.m_grids[PlasmaDomain::d_x];
    const Grid& d_y = m_pd.m_grids[PlasmaDomain::d_y];
    const Grid& be_x = m_pd.m_grids[PlasmaDomain::be_x];
    const Grid& be_y = m_pd.m_grids[PlasmaDomain::be_y];
    const Grid& be_z = m_pd.m_grids[PlasmaDomain::be_z];
    // continuity equations
    std::vector<Grid> v = {grids[v_x],grids[v_y]};
    Grid d_rho_dt = -m_pd.transportDivergence2D(grids[rho],v);
    // magnetic forces
    Grid curl_db = m_pd.curl2D(grids[bi_x],grids[bi_y])/(4.0*PI);
    std::vector<Grid> external_mag_force = Grid::CrossProductZ2D(curl_db,{be_x,be_y});
    std::vector<Grid> internal_mag_force = Grid::CrossProductZ2D(curl_db,{grids[bi_x],grids[bi_y]});
    std::vector<Grid> internal_z_mag_force = Grid::CrossProduct2DZ(m_pd.curlZ(grids[bi_z]),grids[bi_z]/(4.0*PI)); //force in xy-plane due to B_i field in z direction
    std::vector<Grid> external_z_mag_force = Grid::CrossProduct2DZ(m_pd.curlZ(grids[bi_z]),be_z/(4.0*PI)); //force in xy-plane due to B_e field in z direction
    Grid mag_force_z_external = Grid::CrossProduct2D(m_pd.curlZ(grids[bi_z]),{be_x, be_y})/(4.0*PI); //force in z direction due to B field
    Grid mag_force_z_internal = Grid::CrossProduct2D(m_pd.curlZ(grids[bi_z]),{grids[bi_x], grids[bi_y]})/(4.0*PI); //force in z direction due to B field
    // momentum equations   
    Grid d_mom_x_dt =   - m_pd.transportDivergence2D(grids[mom_x], v)
                        - m_pd.derivative1D(grids[press], 0)
                        + grids[rho]*grids[grav_x]
                        + external_mag_force[0] + internal_mag_force[0]
                        + internal_z_mag_force[0] + external_z_mag_force[0];
    Grid d_mom_y_dt =   - m_pd.transportDivergence2D(grids[mom_y], v)
                        - m_pd.derivative1D(grids[press], 1)
                        + grids[rho]*grids[grav_y]
                        + external_mag_force[1] + internal_mag_force[1]
                        + internal_z_mag_force[1] + external_z_mag_force[1];
    Grid d_mom_z_dt =   - m_pd.transportDivergence2D(grids[mom_z], v)
                        + mag_force_z_external + mag_force_z_internal;
    // energy equations
    Grid d_thermal_energy_dt =  - m_pd.transportDivergence2D(grids[thermal_energy],v)
                                - grids[press]*m_pd.divergence2D(v);
    // induction equations
    // std::vector<Grid> induction_rhs_external = m_pd.curlZ(Grid::CrossProduct2D(v,{be_x,be_y}));
    // std::vector<Grid> induction_rhs_internal = m_pd.curlZ(Grid::CrossProduct2D(v,{grids[bi_x],grids[bi_y]}));
    // std::vector<Grid> vz_cross_be = Grid::CrossProductZ2D(grids[v_z],{be_x,be_y});
    // std::vector<Grid> vz_cross_bi = Grid::CrossProductZ2D(grids[v_z],{grids[bi_x],grids[bi_y]});
    // std::vector<Grid> vxy_cross_bi = Grid::CrossProduct2DZ(v,grids[bi_z]);
    // Grid induction_rhs_z_external = m_pd.curl2D(vz_cross_be[0],vz_cross_be[1]);
    // Grid induction_rhs_z_internal = m_pd.curl2D(vz_cross_bi[0],vz_cross_bi[1]);
    // Grid induction_rhs_z_z = m_pd.curl2D(vxy_cross_bi[0],vxy_cross_bi[1]);
    // Grid d_bi_x_dt = induction_rhs_external[0] + induction_rhs_internal[0];
    // Grid d_bi_y_dt = induction_rhs_external[1] + induction_rhs_internal[1];
    // Grid d_bi_z_dt = induction_rhs_z_external + induction_rhs_z_internal + induction_rhs_z_z;
    
    Grid d_bi_x_dt = - m_pd.transportDivergence2D(grids[bi_x],v) - m_pd.transportDivergence2D(be_x,v)
                     + (grids[bi_x] + be_x)*m_pd.derivative1D(grids[v_x],0)
                     + (grids[bi_y] + be_y)*m_pd.derivative1D(grids[v_x],1);
    Grid d_bi_y_dt = - m_pd.transportDivergence2D(grids[bi_y],v) - m_pd.transportDivergence2D(be_y,v)
                     + (grids[bi_x] + be_x)*m_pd.derivative1D(grids[v_y],0)
                     + (grids[bi_y] + be_y)*m_pd.derivative1D(grids[v_y],1);
    Grid d_bi_z_dt = - m_pd.transportDivergence2D(grids[bi_z],v) - m_pd.transportDivergence2D(be_z,v)
                     + (grids[bi_x] + be_x)*m_pd.derivative1D(grids[v_z],0)
                     + (grids[bi_y] + be_y)*m_pd.derivative1D(grids[v_z],1);
    
    // characteristic boundary evolution
    double global_visc_coeff = m_global_viscosity*0.5*((d_x.square()+d_y.square())/grids[dt]).min(m_pd.m_xl,m_pd.m_yl,m_pd.m_xu,m_pd.m_yu);
    std::vector<Grid> char_evolution = computeTimeDerivativesCharacteristicBoundary(grids,
                                                m_pd.x_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC,
                                                m_pd.x_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC,
                                                m_pd.y_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC,
                                                m_pd.y_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC, global_visc_coeff);

    // compute final time derivatives

    std::vector<Grid> grids_dt {d_rho_dt,d_mom_x_dt,d_mom_y_dt,d_mom_z_dt,d_thermal_energy_dt,d_bi_x_dt,d_bi_y_dt,d_bi_z_dt};
    for (int i=0; i<evolved_variables().size(); i++){
        grids_dt[i] *= m_pd.m_ghost_zone_mask;
        grids_dt[i] += char_evolution[i];
    }
    return grids_dt;
}

// void IdealMHD::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){
//     assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
//     grids[rho] += step*time_derivatives[0];
//     grids[mom_x] += step*time_derivatives[1];
//     grids[mom_y] += step*time_derivatives[2];
//     grids[mom_z] += step*time_derivatives[3];
//     grids[thermal_energy] += step*time_derivatives[4];
//     grids[bi_x] += step*time_derivatives[5];
//     grids[bi_y] += step*time_derivatives[6];
//     grids[bi_z] += step*time_derivatives[7];
//     if(moc_mom_limiting) applyMomThresholdingMoC(grids);
//     if(moc_b_limiting) applyBThresholdingMoC(grids);
//     propagateChanges(grids);
// }


void IdealMHD::applyBThresholdingMoC(std::vector<Grid> &grids){
    std::vector<int> be = {PlasmaDomain::be_x,PlasmaDomain::be_y,PlasmaDomain::be_z};
    std::vector<int> bi = {bi_x,bi_y,bi_z};
    for(int n : {0,1,2}){
        Grid& this_be = m_pd.m_grids[be[n]];
        Grid& this_bi = grids[bi[n]];
        if(m_pd.x_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int j=0; j<m_pd.m_ydim; j++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(this_be(N_GHOST+1,j)+this_bi(N_GHOST+1,j) >= 0.0) this_bi(k,j) = std::max(this_bi(k,j),moc_b_lim_threshold*(this_be(N_GHOST+1,j)+this_bi(N_GHOST+1,j)) - this_be(k,j));
                    else this_bi(k,j) = std::min(this_bi(k,j),moc_b_lim_threshold*(this_be(N_GHOST+1,j)+this_bi(N_GHOST+1,j)) - this_be(k,j));
                }
            }
        }
        if(m_pd.x_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int j=0; j<m_pd.m_ydim; j++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(this_be(m_pd.m_xdim-2-N_GHOST,j)+this_bi(m_pd.m_xdim-2-N_GHOST,j) >= 0.0) this_bi(m_pd.m_xdim-1-k,j) = std::max(this_bi(m_pd.m_xdim-1-k,j),moc_b_lim_threshold*(this_be(m_pd.m_xdim-2-N_GHOST,j)+this_bi(m_pd.m_xdim-2-N_GHOST,j))-this_be(m_pd.m_xdim-1-k,j));
                    else this_bi(m_pd.m_xdim-1-k,j) = std::min(this_bi(m_pd.m_xdim-1-k,j),moc_b_lim_threshold*(this_be(m_pd.m_xdim-2-N_GHOST,j)+this_bi(m_pd.m_xdim-2-N_GHOST,j))-this_be(m_pd.m_xdim-1-k,j));
                }
            }
        }
        if(m_pd.y_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int i=0; i<m_pd.m_xdim; i++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(this_be(i,N_GHOST+1)+this_bi(i,N_GHOST+1) >= 0.0) this_bi(i,k) = std::max(this_bi(i,k),moc_b_lim_threshold*(this_be(i,N_GHOST+1)+this_bi(i,N_GHOST+1))-this_be(i,k));
                    else this_bi(i,k) = std::min(this_bi(i,k),moc_b_lim_threshold*(this_be(i,N_GHOST+1)+this_bi(i,N_GHOST+1))-this_be(i,k));
                }
            }
        }
        if(m_pd.y_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int i=0; i<m_pd.m_xdim; i++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(this_be(i,m_pd.m_xdim-2-N_GHOST)+this_bi(i,m_pd.m_xdim-2-N_GHOST) >= 0.0) this_bi(i,m_pd.m_ydim-1-k) = std::max(this_bi(i,m_pd.m_ydim-1-k),moc_b_lim_threshold*(this_be(i,m_pd.m_xdim-2-N_GHOST)+this_bi(i,m_pd.m_xdim-2-N_GHOST))-this_be(i,m_pd.m_ydim-1-k));
                    else this_bi(i,m_pd.m_ydim-1-k) = std::min(this_bi(i,m_pd.m_ydim-1-k),moc_b_lim_threshold*(this_be(i,m_pd.m_xdim-2-N_GHOST)+this_bi(i,m_pd.m_xdim-2-N_GHOST))-this_be(i,m_pd.m_ydim-1-k));
                }
            }
        }
    }
}

void IdealMHD::applyMomThresholdingMoC(std::vector<Grid> &grids){
    for(int v : {mom_x,mom_y,mom_z}){
        if(m_pd.x_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int j=0; j<m_pd.m_ydim; j++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(grids[v](N_GHOST+1,j) >= 0.0) grids[v](k,j) = std::max(grids[v](k,j),moc_mom_lim_threshold*grids[v](N_GHOST+1,j));
                    else grids[v](k,j) = std::min(grids[v](k,j),moc_mom_lim_threshold*grids[v](N_GHOST+1,j));
                }
            }
        }
        if(m_pd.x_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int j=0; j<m_pd.m_ydim; j++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(grids[v](m_pd.m_xdim-2-N_GHOST,j) >= 0.0) grids[v](m_pd.m_xdim-1-k,j) = std::max(grids[v](m_pd.m_xdim-1-k,j),moc_mom_lim_threshold*grids[v](m_pd.m_xdim-2-N_GHOST,j));
                    else grids[v](m_pd.m_xdim-1-k,j) = std::min(grids[v](m_pd.m_xdim-1-k,j),moc_mom_lim_threshold*grids[v](m_pd.m_xdim-2-N_GHOST,j));
                }
            }
        }
        if(m_pd.y_bound_1==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int i=0; i<m_pd.m_xdim; i++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(grids[v](i,N_GHOST+1) >= 0.0) grids[v](i,k) = std::max(grids[v](i,k),moc_mom_lim_threshold*grids[v](i,N_GHOST+1));
                    else grids[v](i,k) = std::min(grids[v](i,k),moc_mom_lim_threshold*grids[v](i,N_GHOST+1));
                }
            }
        }
        if(m_pd.y_bound_2==PlasmaDomain::BoundaryCondition::OpenMoC){
            for(int i=0; i<m_pd.m_xdim; i++){
                for(int k=0; k<N_GHOST+1; k++){
                    if(grids[v](i,m_pd.m_ydim-2-N_GHOST) >= 0.0) grids[v](i,m_pd.m_ydim-1-k) = std::max(grids[v](i,m_pd.m_ydim-1-k),moc_mom_lim_threshold*grids[v](i,m_pd.m_ydim-2-N_GHOST));
                    else grids[v](i,m_pd.m_ydim-1-k) = std::min(grids[v](i,m_pd.m_ydim-1-k),moc_mom_lim_threshold*grids[v](i,m_pd.m_ydim-2-N_GHOST));
                }
            }
        }
    }
}


void IdealMHD::recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) const
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[press] = 2*grids[n]*K_B*grids[temp].max(m_pd.temp_min);
    grids[thermal_energy] = grids[press]/(m_pd.m_adiabatic_index - 1.0);
}

void IdealMHD::recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) const
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[rho] = grids[n]*m_pd.m_ion_mass;
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[v_z] = grids[mom_z]/grids[rho];
    grids[kinetic_energy] = 0.5*grids[rho]*(grids[v_x].square()+grids[v_y].square());
    grids[thermal_energy] = grids[thermal_energy].max(m_pd.thermal_energy_min);
    grids[press] = (m_pd.m_adiabatic_index - 1.0)*grids[thermal_energy];
    grids[temp] = (grids[press]/(2*K_B*grids[n])).max(m_pd.temp_min);
    grids[b_x] = m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x];
    grids[b_y] = m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_z] = m_pd.m_grids[PlasmaDomain::be_z] + grids[bi_z];
    grids[b_mag] = (grids[b_x].square() + grids[b_y].square() + grids[b_z].square()).sqrt();
    grids[b_hat_x] = grids[b_x]/grids[b_mag];
    grids[b_hat_y] = grids[b_y]/grids[b_mag];
    grids[b_hat_z] = grids[b_z]/grids[b_mag];
    catchNullFieldDirection(grids);
}

void IdealMHD::catchNullFieldDirection(std::vector<Grid> &grids) const
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(grids[b_mag](i,j) == 0.0){
                grids[b_hat_x](i,j) = 0.0;
                grids[b_hat_y](i,j) = 0.0;
                grids[b_hat_z](i,j) = 0.0;
            }
        }
    }
}

void IdealMHD::recomputeDT(std::vector<Grid>& grids) const
{
    Grid c_s = (m_pd.m_adiabatic_index*grids[press]/grids[rho]).sqrt();
    Grid c_s_sq = c_s.square();

    Grid v_alfven = grids[b_mag]/(4.0*PI*grids[rho]).sqrt();
    Grid v_alfven_sq = v_alfven.square();
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);

    //Magnetoacoustic modes
    Grid delta = (one - 4.0*c_s_sq*v_alfven_sq/(c_s_sq+v_alfven_sq).square()).sqrt();
    Grid v_fast = (0.5*(c_s_sq + v_alfven_sq)*(one + delta)).sqrt();
    Grid v_slow = (0.5*(c_s_sq + v_alfven_sq)*(one - delta)).sqrt();

    //Bulk velocity transit time
    Grid diagonals = (m_pd.m_grids[PlasmaDomain::d_x].square() + m_pd.m_grids[PlasmaDomain::d_y].square()).sqrt();
    Grid vel_mag = (grids[v_x].square() + grids[v_y].square()).sqrt();

    // m_grids[dt] = diagonals/(c_s + v_alfven + v_fast + v_slow + vel_mag);
    grids[dt] = diagonals/(vel_mag + c_s.max(v_alfven).max(v_fast).max(v_slow));
}

// Returns time derivative from characteristic boundary cond for the quantities
// {rho, mom_x, mom_y, thermal_energy, b_x, b_y} in that order
std::vector<Grid> IdealMHD::computeTimeDerivativesCharacteristicBoundary(const std::vector<Grid> &grids, bool x_bound_1, bool x_bound_2, bool y_bound_1, bool y_bound_2, double visc_coeff) const
{
    std::vector<std::vector<Grid> > x_lower, x_upper, y_lower, y_upper;
    if(x_bound_1) x_lower = singleBoundaryTermsMOC(grids, 0, true, visc_coeff);
    if(x_bound_2) x_upper = singleBoundaryTermsMOC(grids, 0, false, visc_coeff);
    if(y_bound_1) y_lower = singleBoundaryTermsMOC(grids, 1, true, visc_coeff);
    if(y_bound_2) y_upper = singleBoundaryTermsMOC(grids, 1, false, visc_coeff);
    int num_vars = 8;
    std::vector<bool> bound_flags = {x_bound_1,x_bound_2,y_bound_1,y_bound_2};
    std::vector<std::vector<std::vector<Grid> > > bound_vals = {x_lower, x_upper, y_lower, y_upper};

    std::vector<Grid> results(num_vars,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    for(int i=0; i<num_vars; i++){
        for(int k=0; k<bound_flags.size(); k++) if(bound_flags[k]) results[i] += bound_vals[k][i][0] + bound_vals[k][i][1];
    }

    Grid mom_x_term = grids[rho]*results[1] + grids[v_x]*results[0];
    Grid mom_y_term = grids[rho]*results[2] + grids[v_y]*results[0];
    Grid mom_z_term = grids[rho]*results[3] + grids[v_z]*results[0];
    results[1] = mom_x_term;
    results[2] = mom_y_term;
    results[3] = mom_z_term;

    return results;
}

// Returns {normal term, parallel term} (in full-sized grid) for single MOC boundary,
// corresponding to each of {rho, v_x, v_y, b_x, b_y, thermal_energy}
std::vector<std::vector<Grid>> IdealMHD::singleBoundaryTermsMOC(const std::vector<Grid> &grids, int boundary_index, bool boundary_lower, double visc_coeff) const
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    std::vector<int> x_idx, y_idx; //Full extent of boundary's ghost zone region
    std::vector<int> x_idx_interior, y_idx_interior; //Full extent of boundary's ghost zone not overlapping with other ghost zones
    std::vector<int> x_idx_padded, y_idx_padded; //Full extent of boundary's ghost zone region, reduced for central differencing safety
    std::vector<int> x_idx_lower_corner, x_idx_upper_corner, y_idx_lower_corner, y_idx_upper_corner;
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
    y_idx_interior = y_idx;
    y_idx_padded = y_idx;
    x_idx_lower_corner = x_idx;
    x_idx_upper_corner = x_idx;
    y_idx_lower_corner = y_idx;
    y_idx_upper_corner = y_idx;
    if(boundary_index == 0 && m_pd.y_bound_1 != PlasmaDomain::BoundaryCondition::Periodic){
        y_idx_interior = {y_idx[0]+N_GHOST,y_idx.back()-N_GHOST};
        y_idx_padded = {y_idx[0]+1,y_idx.back()-1};
        if(m_pd.y_bound_1 != PlasmaDomain::BoundaryCondition::OpenMoC){
            y_idx[0] = y_idx_interior[0];
            y_idx_padded[0] = y_idx_interior[0];
        }
        if(m_pd.y_bound_2 != PlasmaDomain::BoundaryCondition::OpenMoC){
            y_idx[1] = y_idx_interior[1];
            y_idx_padded[1] = y_idx_interior[1];
        }
        y_idx_lower_corner = {y_idx[0],y_idx_padded[0]-1};
        y_idx_upper_corner = {y_idx_padded[1]+1,y_idx[1]};
    } else if(boundary_index == 1 && m_pd.x_bound_1 != PlasmaDomain::BoundaryCondition::Periodic){
        x_idx_interior = {x_idx[0]+N_GHOST,x_idx.back()-N_GHOST};
        x_idx_padded = {x_idx[0]+1,x_idx.back()-1};
        if(m_pd.x_bound_1 != PlasmaDomain::BoundaryCondition::OpenMoC){
            x_idx[0] = x_idx_interior[0];
            x_idx_padded[0] = x_idx_interior[0];
        }
        if(m_pd.x_bound_2 != PlasmaDomain::BoundaryCondition::OpenMoC){
            x_idx[1] = x_idx_interior[1];
            x_idx_padded[1] = x_idx_interior[1];
        }
        x_idx_lower_corner = {x_idx[0],x_idx_padded[0]-1};
        x_idx_upper_corner = {x_idx_padded[1]+1,x_idx[1]};
    }

    // Mask grid where equal to 1 inside of boundary ghost zone, 0 elsewhere
    Grid boundary_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid interior_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid lower_corner_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid upper_corner_mask = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    for(int i=x_idx[0]; i<=x_idx[1]; i++) for(int j=y_idx[0]; j<=y_idx[1]; j++) boundary_mask(i,j) = 1.0;
    for(int i=x_idx_interior[0]; i<=x_idx_interior[1]; i++) for(int j=y_idx_interior[0]; j<=y_idx_interior[1]; j++) interior_mask(i,j) = 1.0;

    std::vector<int> vels = {v_x,v_y,v_z};
    std::vector<int> bes = {PlasmaDomain::be_x, PlasmaDomain::be_y, PlasmaDomain::be_z};
    std::vector<int> bis = {bi_x, bi_y, bi_z};

    int v_perp = vels[boundary_index];
    int v_para = vels[(boundary_index + 1)%2];
    int v_guide = v_z;
    int be_perp = bes[boundary_index];
    int be_para = bes[(boundary_index + 1)%2];
    int be_guide = PlasmaDomain::be_z;
    int bi_perp = bis[boundary_index];
    int bi_para = bis[(boundary_index + 1)%2];
    int bi_guide = bi_z;
    int grav_perp = (boundary_index == 0) ? grav_x : grav_y;
    int grav_para = (boundary_index == 1) ? grav_x : grav_y;
    bool positive_forward = !boundary_lower;

    const Grid& m_be_perp = m_pd.m_grids[be_perp];
    const Grid& m_be_para = m_pd.m_grids[be_para];
    const Grid& m_be_guide = m_pd.m_grids[be_guide];
    const Grid& m_bi_perp = grids[bi_perp];
    const Grid& m_bi_para = grids[bi_para];
    const Grid& m_bi_guide = grids[bi_guide];

    Grid m_b_perp = m_be_perp + m_bi_perp;
    Grid m_b_para = m_be_para + m_bi_para;
    Grid m_b_guide = m_be_guide + m_bi_guide;

    Grid m_b_h_mag = (m_b_para.square() + m_b_guide.square()).sqrt();

    // double b_guide = m_guide_field; //constant guide field in third dimension
    Grid s_perp, R_para, R_guide, c_s_sq, c_a_sq, c_perp_sq, c_plus_sq, c_minus_sq, alpha_plus_sq, alpha_minus_sq, c_s, c_perp, c_plus, c_minus, alpha_plus, alpha_minus;
    s_perp=R_para=R_guide=c_a_sq=c_perp_sq=c_plus_sq=c_minus_sq=alpha_plus_sq=alpha_minus_sq=c_perp=c_plus=c_minus=alpha_plus=alpha_minus = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    c_s = c_s_sq = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim); //set to one to avoid problems when dividing by c_s_sq
    for(int i = x_idx[0]; i <= x_idx[1]; i++){
        for(int j = y_idx[0]; j <= y_idx[1]; j++){
            s_perp(i,j) = m_b_perp(i,j) > 0.0 ? 1.0 : -1.0; ///std::abs(m_b_perp(i,j));
            R_para(i,j) = m_b_para(i,j)/m_b_h_mag(i,j);
            R_guide(i,j) = m_b_guide(i,j)/m_b_h_mag(i,j);
            c_s_sq(i,j) = m_pd.m_adiabatic_index*grids[press](i,j)/grids[rho](i,j);
            c_a_sq(i,j) = grids[b_mag](i,j)*grids[b_mag](i,j)/(4.0*PI*grids[rho](i,j));
            c_perp_sq(i,j) = m_b_perp(i,j)*m_b_perp(i,j)/(4.0*PI*grids[rho](i,j));
            c_plus_sq(i,j) = (c_a_sq(i,j) + c_s_sq(i,j))/2.0 + std::sqrt(std::abs((std::pow(c_a_sq(i,j) + c_s_sq(i,j),2.0)/4.0 - c_perp_sq(i,j)*c_s_sq(i,j))));
            c_minus_sq(i,j) = (c_a_sq(i,j) + c_s_sq(i,j))/2.0 - std::sqrt(std::abs((std::pow(c_a_sq(i,j) + c_s_sq(i,j),2.0)/4.0 - c_perp_sq(i,j)*c_s_sq(i,j))));
            alpha_plus_sq(i,j) = std::max((c_s_sq(i,j) - c_minus_sq(i,j))/(c_plus_sq(i,j) - c_minus_sq(i,j)),0.0);
            alpha_minus_sq(i,j) = std::max((c_plus_sq(i,j) - c_s_sq(i,j))/(c_plus_sq(i,j) - c_minus_sq(i,j)),0.0);
            c_s(i,j) = std::sqrt(c_s_sq(i,j));
            c_perp(i,j) = std::sqrt(c_perp_sq(i,j));
            c_plus(i,j) = std::sqrt(c_plus_sq(i,j)), c_minus(i,j) = std::sqrt(c_minus_sq(i,j));
            alpha_plus(i,j) = std::sqrt(alpha_plus_sq(i,j)), alpha_minus(i,j) = std::sqrt(alpha_minus_sq(i,j));
        }
    }

    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(m_b_perp(i,j) == 0.0) s_perp(i,j) = 0.0;
        }
    }

    //Compute ordinary perpendicular terms
    Grid b_para_grad = m_pd.characteristicBartonDerivative1D(m_be_para,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
                        + m_pd.characteristicBartonDerivative1D(m_bi_para,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back());
    Grid b_guide_grad = m_pd.characteristicBartonDerivative1D(m_be_guide,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
                        + m_pd.characteristicBartonDerivative1D(m_bi_guide,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back());
    Grid v_para_grad = m_pd.characteristicBartonDerivative1D(grids[v_para],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back());
    Grid v_perp_grad = m_pd.characteristicBartonDerivative1D(grids[v_perp],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back());
    Grid v_guide_grad = m_pd.characteristicBartonDerivative1D(grids[v_guide],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back());
    Grid press_grad = m_pd.characteristicBartonDerivative1D(grids[press],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back());
    
    Grid d_1 = grids[v_perp]*(
        m_pd.characteristicBartonDerivative1D(m_b_perp,positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );
    Grid d_2 = grids[v_perp]*(
        (grids[thermal_energy] + grids[press])*m_pd.characteristicBartonDerivative1D(grids[rho],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
        - grids[rho]*m_pd.characteristicBartonDerivative1D(grids[thermal_energy],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back())
    );
    Grid d_3 = s_perp*(boundary_index == 0 ? 1.0 : -1.0)*(grids[v_perp] + c_perp)*(
        -s_perp*R_guide*v_para_grad + s_perp*R_para*v_guide_grad
        + R_guide/(4.0*PI*grids[rho]).sqrt()*b_para_grad - R_para/(4.0*PI*grids[rho]).sqrt()*b_guide_grad
    );
    Grid d_4 = s_perp*(boundary_index == 0 ? 1.0 : -1.0)*(grids[v_perp] - c_perp)*(
        s_perp*R_guide*v_para_grad - s_perp*R_para*v_guide_grad
        + R_guide/(4.0*PI*grids[rho]).sqrt()*b_para_grad - R_para/(4.0*PI*grids[rho]).sqrt()*b_guide_grad
    );
    Grid d_5 = (grids[v_perp] + c_plus)*(
        alpha_plus/grids[rho]*press_grad - s_perp*R_para*c_minus*alpha_minus*v_para_grad - s_perp*R_guide*c_minus*alpha_minus*v_guide_grad
        + c_plus*alpha_plus*v_perp_grad + R_para*c_s*alpha_minus/(4.0*PI*(grids[rho])).sqrt()*b_para_grad + R_guide*c_s*alpha_minus/(4.0*PI*(grids[rho])).sqrt()*b_guide_grad
    );
    Grid d_6 = (grids[v_perp] - c_plus)*(
        alpha_plus/grids[rho]*press_grad + s_perp*R_para*c_minus*alpha_minus*v_para_grad + s_perp*R_guide*c_minus*alpha_minus*v_guide_grad
        - c_plus*alpha_plus*v_perp_grad + R_para*c_s*alpha_minus/(4.0*PI*(grids[rho])).sqrt()*b_para_grad + R_guide*c_s*alpha_minus/(4.0*PI*(grids[rho])).sqrt()*b_guide_grad
    );
    Grid d_7 = (grids[v_perp] + c_minus)*(
        alpha_minus/grids[rho]*press_grad + s_perp*R_para*c_plus*alpha_plus*v_para_grad + s_perp*R_guide*c_plus*alpha_plus*v_guide_grad
        + c_minus*alpha_minus*v_perp_grad - R_para*c_s*alpha_plus/(4.0*PI*(grids[rho])).sqrt()*b_para_grad - R_guide*c_s*alpha_plus/(4.0*PI*(grids[rho])).sqrt()*b_guide_grad
    );
    Grid d_8 = (grids[v_perp] - c_minus)*(
        alpha_minus/grids[rho]*press_grad - s_perp*R_para*c_plus*alpha_plus*v_para_grad - s_perp*R_guide*c_plus*alpha_plus*v_guide_grad
        - c_minus*alpha_minus*v_perp_grad - R_para*c_s*alpha_plus/(4.0*PI*(grids[rho])).sqrt()*b_para_grad - R_guide*c_s*alpha_plus/(4.0*PI*(grids[rho])).sqrt()*b_guide_grad
    );
    int para_index = boundary_index == 0 ? 1 : 0;

    //Compute alternate perpendicular terms (for negating inflowing characteristics)
    Grid b_perp_sq_grad_alt = 2.0*m_b_perp*(
                        m_pd.derivative1D(m_be_perp,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1])
                        +m_pd.derivative1DBackward(m_be_perp,true,para_index,x_idx_upper_corner[0],y_idx_upper_corner[0],x_idx_upper_corner[1],y_idx_upper_corner[1])
                        +m_pd.derivative1DBackward(m_be_perp,false,para_index,x_idx_lower_corner[0],y_idx_lower_corner[0],x_idx_lower_corner[1],y_idx_lower_corner[1])
                        +m_pd.derivative1D(m_bi_perp,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1])
                        +m_pd.derivative1DBackward(m_bi_perp,true,para_index,x_idx_upper_corner[0],y_idx_upper_corner[0],x_idx_upper_corner[1],y_idx_upper_corner[1])
                        +m_pd.derivative1DBackward(m_bi_perp,false,para_index,x_idx_lower_corner[0],y_idx_lower_corner[0],x_idx_lower_corner[1],y_idx_lower_corner[1]));
    Grid press_grad_alt = m_pd.derivative1D(grids[press],para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1])
                +m_pd.derivative1DBackward(grids[press],true,para_index,x_idx_upper_corner[0],y_idx_upper_corner[0],x_idx_upper_corner[1],y_idx_upper_corner[1])
                +m_pd.derivative1DBackward(grids[press],false,para_index,x_idx_lower_corner[0],y_idx_lower_corner[0],x_idx_lower_corner[1],y_idx_lower_corner[1]);
    Grid b_perp_grad_alt = m_pd.derivative1D(m_be_perp,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1])
            +m_pd.derivative1DBackward(m_be_perp,true,para_index,x_idx_upper_corner[0],y_idx_upper_corner[0],x_idx_upper_corner[1],y_idx_upper_corner[1])
            +m_pd.derivative1DBackward(m_be_perp,false,para_index,x_idx_lower_corner[0],y_idx_lower_corner[0],x_idx_lower_corner[1],y_idx_lower_corner[1])
            +m_pd.derivative1D(m_bi_perp,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1])
            +m_pd.derivative1DBackward(m_bi_perp,true,para_index,x_idx_upper_corner[0],y_idx_upper_corner[0],x_idx_upper_corner[1],y_idx_upper_corner[1])
            +m_pd.derivative1DBackward(m_bi_perp,false,para_index,x_idx_lower_corner[0],y_idx_lower_corner[0],x_idx_lower_corner[1],y_idx_lower_corner[1]);
    Grid b_guide_grad_alt = m_pd.derivative1D(m_be_guide,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1])
                    +m_pd.derivative1DBackward(m_be_guide,true,para_index,x_idx_upper_corner[0],y_idx_upper_corner[0],x_idx_upper_corner[1],y_idx_upper_corner[1])
                    +m_pd.derivative1DBackward(m_be_guide,false,para_index,x_idx_lower_corner[0],y_idx_lower_corner[0],x_idx_lower_corner[1],y_idx_lower_corner[1])
                    +m_pd.derivative1D(m_bi_guide,para_index,x_idx_padded[0],y_idx_padded[0],x_idx_padded[1],y_idx_padded[1])
                    +m_pd.derivative1DBackward(m_bi_guide,true,para_index,x_idx_upper_corner[0],y_idx_upper_corner[0],x_idx_upper_corner[1],y_idx_upper_corner[1])
                    +m_pd.derivative1DBackward(m_bi_guide,false,para_index,x_idx_lower_corner[0],y_idx_lower_corner[0],x_idx_lower_corner[1],y_idx_lower_corner[1]);

    Grid d_3_alt = s_perp*s_perp*(boundary_index == 0 ? 1.0 : -1.0)*( (-1.0*interior_mask*R_guide*grids[grav_para]
        + (R_guide*(press_grad_alt + b_perp_sq_grad_alt/(8.0*PI)) + m_b_h_mag*(b_guide_grad_alt)/(4.0*PI))/grids[rho]) ); //R_para*m_b_para*(b_guide_grad_alt)/(4.0*PI))/grids[rho]) );
    Grid d_5_alt = c_plus*alpha_plus*(interior_mask*grids[grav_perp] + m_b_para*b_perp_grad_alt/(4.0*PI*grids[rho]))
        + c_minus*alpha_minus*s_perp*(-interior_mask*R_para*grids[grav_para] + R_para*(press_grad_alt+b_perp_sq_grad_alt/(8.0*PI))/grids[rho]);
    Grid d_7_alt = c_minus*alpha_minus*(interior_mask*grids[grav_perp] + m_b_para*b_perp_grad_alt/(4.0*PI*grids[rho]))
        - c_plus*alpha_plus*s_perp*(-interior_mask*R_para*grids[grav_para] + R_para*(press_grad_alt+b_perp_sq_grad_alt/(8.0*PI))/grids[rho]);

    //Check for inflowing characteristic speeds and swap out alternate values as necessary
    int inflow_dir = boundary_lower ? 1 : -1;
    for(int i = x_idx[0]; i <= x_idx[1]; i++){
        for(int j = y_idx[0]; j <= y_idx[1]; j++){
            if(inflow_dir*grids[v_perp](i,j) > 0.0) { d_1(i,j) = 0.0; d_2(i,j) = 0.0;}
            if(inflow_dir*(grids[v_perp](i,j) + c_perp(i,j)) > 0.0) {d_3(i,j) = d_3_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) - c_perp(i,j)) > 0.0) {d_4(i,j) = -d_3_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) + c_plus(i,j)) > 0.0) {d_5(i,j) = d_5_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) - c_plus(i,j)) > 0.0) {d_6(i,j) = -d_5_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) + c_minus(i,j)) > 0.0) {d_7(i,j) = d_7_alt(i,j);}
            if(inflow_dir*(grids[v_perp](i,j) - c_minus(i,j)) > 0.0) {d_8(i,j) = -d_7_alt(i,j);}
        }
    }
    //{rho, v_x, v_y, b_x, b_y, thermal_energy}
    std::vector<Grid> rho_terms, v_perp_terms, v_para_terms, v_guide_terms, b_perp_terms, b_para_terms, b_guide_terms, thermal_energy_terms;

    // Compute all normal terms
    rho_terms.push_back(-((m_pd.m_adiabatic_index/grids[rho])*d_2 + 0.5*grids[rho]*alpha_plus*(d_5+d_6) + 0.5*grids[rho]*alpha_minus*(d_7+d_8))/c_s_sq);
    v_para_terms.push_back(-0.5*s_perp*((boundary_index == 0 ? 1.0 : -1.0)*R_guide*(-d_3 + d_4)
                            + (c_minus*alpha_minus/c_s_sq*R_para*(-d_5 + d_6) + c_plus*alpha_plus/c_s_sq*R_para*(d_7 - d_8)))//);
                            + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_para],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back()));
    v_guide_terms.push_back(-0.5*s_perp*((boundary_index == 0 ? 1.0 : -1.0)*R_para*(d_3 - d_4)
                            + (c_minus*alpha_minus/c_s_sq*R_guide*(-d_5 + d_6) + c_plus*alpha_plus/c_s_sq*R_guide*(d_7 - d_8)))//);
                            + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_guide],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back()));
    v_perp_terms.push_back(-0.5/c_s_sq*(c_plus*alpha_plus*(d_5 - d_6) + c_minus*alpha_minus*(d_7 - d_8))//);
                            + visc_coeff*m_pd.secondDerivative1DBackward(grids[v_perp],positive_forward,boundary_index,x_idx[0],y_idx[0],x_idx.back(),y_idx.back()));
    thermal_energy_terms.push_back(-(0.5*(grids[thermal_energy]+grids[press])*alpha_plus*(d_5+d_6) + 0.5*(grids[thermal_energy]+grids[press])*alpha_minus*(d_7+d_8))/c_s_sq);
    b_para_terms.push_back(-(PI*grids[rho]).sqrt()*( (boundary_index == 0 ? 1.0 : -1.0)*R_guide*(d_3 + d_4) + alpha_minus/c_s*R_para*(d_5 + d_6) - alpha_plus/c_s*R_para*(d_7 + d_8) ));
    b_guide_terms.push_back(-(PI*grids[rho]).sqrt()*( (boundary_index == 0 ? 1.0 : -1.0)*-1.0*R_para*(d_3 + d_4) + alpha_minus/c_s*R_guide*(d_5 + d_6) - alpha_plus/c_s*R_guide*(d_7 + d_8) ));
    b_perp_terms.push_back(-d_1);

    // Compute all parallel terms
    rho_terms.push_back(-m_pd.transportDerivative1D(grids[rho], grids[v_para], para_index,
                        x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    v_para_terms.push_back(-1.0/grids[rho]*(m_pd.derivative1D(grids[press],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                            +2.0*m_b_para/(8.0*PI)*(
                                                m_pd.derivative1D(m_be_para,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                                +m_pd.derivative1D(m_bi_para,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                            )
                                            +2.0*m_b_perp/(8.0*PI)*(
                                                m_pd.derivative1D(m_be_perp,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                                +m_pd.derivative1D(m_bi_perp,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                            )
                                            +2.0*m_b_guide/(8.0*PI)*(
                                                m_pd.derivative1D(m_be_guide,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                                +m_pd.derivative1D(m_bi_guide,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                            )
                                        )
                        -( m_pd.transportDerivative1D(grids[rho]*grids[v_para],grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                            - grids[v_para]*m_pd.transportDerivative1D(grids[rho],grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]) )/grids[rho]
                        + m_b_para/(4.0*PI*grids[rho])*(
                                m_pd.derivative1D(m_be_para,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                +m_pd.derivative1D(m_bi_para,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]) )
                        + interior_mask*grids[grav_para]
                        + visc_coeff*m_pd.secondDerivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    v_guide_terms.push_back( -( m_pd.transportDerivative1D(grids[rho]*grids[v_guide],grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                            - grids[v_guide]*m_pd.transportDerivative1D(grids[rho],grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]) )/grids[rho]
                        +m_b_para/(4.0*PI*grids[rho])*(
                                m_pd.derivative1D(m_be_guide,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                +m_pd.derivative1D(m_bi_guide,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]) )
                        + visc_coeff*m_pd.secondDerivative1D(grids[v_guide],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    v_perp_terms.push_back( -( m_pd.transportDerivative1D(grids[rho]*grids[v_perp],grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                            - grids[v_perp]*m_pd.transportDerivative1D(grids[rho],grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]) )/grids[rho]
                        +m_b_para/(4.0*PI*grids[rho])*(
                                m_pd.derivative1D(m_be_perp,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                +m_pd.derivative1D(m_bi_perp,para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]) )
                        +interior_mask*grids[grav_perp]
                        + visc_coeff*m_pd.secondDerivative1D(grids[v_perp],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    thermal_energy_terms.push_back(-m_pd.transportDerivative1D(grids[thermal_energy],grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                                -grids[press]*m_pd.derivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    b_para_terms.push_back(-(m_pd.transportDerivative1D(m_be_para,grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                              + m_pd.transportDerivative1D(m_bi_para,grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]))
                            +m_b_para*m_pd.derivative1D(grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    b_guide_terms.push_back(-(m_pd.transportDerivative1D(m_be_guide,grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                              + m_pd.transportDerivative1D(m_bi_guide,grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]))
                            +m_b_para*m_pd.derivative1D(grids[v_guide],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));
    b_perp_terms.push_back(-(m_pd.transportDerivative1D(m_be_perp,grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1])
                              + m_pd.transportDerivative1D(m_bi_perp,grids[v_para],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]))
                            +m_b_para*m_pd.derivative1D(grids[v_perp],para_index,x_idx_interior[0],y_idx_interior[0],x_idx_interior[1],y_idx_interior[1]));


    std::vector<Grid> zero(2,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    return {rho_terms,
            v_x == v_para ? v_para_terms : v_perp_terms, v_y == v_para ? v_para_terms : v_perp_terms, v_guide_terms,
            thermal_energy_terms,
            boundary_index == 0 ? b_perp_terms : b_para_terms, boundary_index == 1 ? b_perp_terms : b_para_terms, b_guide_terms};
}

