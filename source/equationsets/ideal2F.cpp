#include "ideal2F.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>

Ideal2F::Ideal2F(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

void Ideal2F::parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for (int i=0; i<lhs.size(); i++){
        if (lhs[i] == "use_sub_cycling") m_use_sub_cycling = (rhs[i] == "true");
        else if (lhs[i] == "epsilon_courant") m_epsilon_courant = stod(rhs[i]);
        else if (lhs[i] == "verbose_2F") m_verbose = (rhs[i] == "true");
        else if (lhs[i] == "viscosity_opt") m_viscosity_opt = rhs[i];
        else if (lhs[i] == "remove_curl_terms") m_remove_curl_terms = (rhs[i] == "true");
        else if (lhs[i] == "global_viscosity") m_global_viscosity = stod(rhs[i]);
        else{
            std::cerr << lhs[i] << " is not recognized for this equation set." << std::endl;
            assert(false);
        }
    }
}

void Ideal2F::setupEquationSetDerived()
{
    // do not subcycle if curl terms are not used
    if (m_remove_curl_terms) m_use_sub_cycling = false;
}

std::vector<Grid> Ideal2F::computeTimeDerivativesDerived(const std::vector<Grid> &grids) const
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // PlasmaDomain grid references for more concise notation
    const Grid& d_x = m_pd.m_grids[PlasmaDomain::d_x];
    const Grid& d_y = m_pd.m_grids[PlasmaDomain::d_y];
    // continuity equations
    std::vector<Grid> v_i = {grids[i_v_x],grids[i_v_y]};
    std::vector<Grid> v_e = {grids[e_v_x],grids[e_v_y]};
    Grid d_i_rho_dt = -m_pd.transportDivergence2D(grids[i_rho],v_i);
    Grid d_e_rho_dt = -m_pd.transportDivergence2D(grids[e_rho],v_e);
    // lorentz forces
    std::vector<Grid> b = {grids[b_x],grids[b_y]};
    std::vector<Grid> i_v_cross_B(2,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    i_v_cross_B[0] = grids[i_v_y]*grids[b_z];
    i_v_cross_B[1] = -grids[i_v_x]*grids[b_z];
    std::vector<Grid> e_v_cross_B(2,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    e_v_cross_B[0] = grids[e_v_y]*grids[b_z];
    e_v_cross_B[1] = -grids[e_v_x]*grids[b_z];
    Grid i_F_x =  E*grids[i_n]*(grids[E_x] + i_v_cross_B[0]/C);
    Grid i_F_y =  E*grids[i_n]*(grids[E_y] + i_v_cross_B[1]/C);
    Grid e_F_x = -E*grids[e_n]*(grids[E_x] + e_v_cross_B[0]/C);
    Grid e_F_y = -E*grids[e_n]*(grids[E_y] + e_v_cross_B[1]/C); 
    // momentum equations   
    Grid d_i_mom_x_dt = - m_pd.transportDivergence2D(grids[i_mom_x], v_i)
                        - m_pd.derivative1D(grids[i_press], 0)
                        + grids[i_rho]*grids[grav_x] + i_F_x;
    Grid d_i_mom_y_dt = - m_pd.transportDivergence2D(grids[i_mom_y], v_i)
                        - m_pd.derivative1D(grids[i_press], 1)
                        + grids[i_rho]*grids[grav_y] + i_F_y;
    Grid d_e_mom_x_dt = - m_pd.transportDivergence2D(grids[e_mom_x], v_e)
                        - m_pd.derivative1D(grids[e_press], 0)
                        + grids[e_rho]*grids[grav_x] + e_F_x;
    Grid d_e_mom_y_dt = - m_pd.transportDivergence2D(grids[e_mom_y], v_e)
                        - m_pd.derivative1D(grids[e_press], 1)
                        + grids[e_rho]*grids[grav_y] + e_F_y;
    // energy equations
    Grid d_i_thermal_energy_dt =   - m_pd.transportDivergence2D(grids[i_thermal_energy],v_i)
                            - grids[i_press]*m_pd.divergence2D(v_i);
    Grid d_e_thermal_energy_dt =   - m_pd.transportDivergence2D(grids[e_thermal_energy],v_e)
                            - grids[e_press]*m_pd.divergence2D(v_e);
    // optionally evolve electromagnetic fields
    Grid d_E_x_dt,d_E_y_dt,d_E_z_dt,d_bi_x_dt,d_bi_y_dt,d_bi_z_dt;
    if (!m_use_sub_cycling){
        if (!m_remove_curl_terms){
            d_E_x_dt = + C*m_pd.derivative1D(grids[b_z],1) - 4.*PI*grids[j_x];
            d_E_y_dt = - C*m_pd.derivative1D(grids[b_z],0) - 4.*PI*grids[j_y];
            d_E_z_dt = + C*(m_pd.derivative1D(grids[b_y],0) - m_pd.derivative1D(grids[b_x],1));
            d_bi_x_dt = - C*m_pd.derivative1D(grids[E_z],1);
            d_bi_y_dt = + C*m_pd.derivative1D(grids[E_z],0);
            d_bi_z_dt = + C*(m_pd.derivative1D(grids[E_x],1) - m_pd.derivative1D(grids[E_y],0));
        }
        else{
            d_E_x_dt = - 4.*PI*grids[j_x];
            d_E_y_dt = - 4.*PI*grids[j_y];
            d_E_z_dt = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
            d_bi_x_dt = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
            d_bi_y_dt = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
            d_bi_z_dt = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
        }
    }

    // apply ghost zone mask
    std::vector<Grid> grids_dt {d_i_rho_dt,d_e_rho_dt,d_i_mom_x_dt,d_i_mom_y_dt,d_e_mom_x_dt,d_e_mom_y_dt,d_i_thermal_energy_dt,d_e_thermal_energy_dt,d_E_x_dt,d_E_y_dt,d_E_z_dt,d_bi_x_dt,d_bi_y_dt,d_bi_z_dt};
    for (int i=0; i<evolved_variables().size(); i++) grids_dt[i] *= m_pd.m_ghost_zone_mask;
    return grids_dt;
}

void Ideal2F::recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) const
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[i_n] = grids[i_rho]/m_pd.m_ion_mass;
    grids[e_n] = grids[e_rho]/M_ELECTRON;
    grids[i_press] = grids[i_n]*K_B*grids[i_temp];
    grids[e_press] = grids[e_n]*K_B*grids[e_temp];
    grids[press] = grids[i_press] + grids[e_press];
    grids[i_thermal_energy] = grids[i_press]/(m_pd.m_adiabatic_index - 1.0);
    grids[e_thermal_energy] = grids[e_press]/(m_pd.m_adiabatic_index - 1.0);
}

void Ideal2F::recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) const
{
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
    grids[n] = grids[i_n];
    grids[dn] = grids[i_n] - grids[e_n];
    grids[divE] = m_pd.divergence2D({grids[E_x],grids[E_y]});
    grids[divB] = m_pd.divergence2D({grids[b_x],grids[b_y]});
    grids[curlE_z] = -(m_pd.derivative1D(grids[E_x],1) - m_pd.derivative1D(grids[E_y],0));
    grids[lapEx] = m_pd.laplacian(grids[E_x]);
    populate_boundary(grids[lapEx]);
    grids[E_ideal] = (M_ELECTRON*m_pd.derivative1D(grids[i_press], 0) - m_pd.m_ion_mass*m_pd.derivative1D(grids[e_press], 0))/(M_ELECTRON+m_pd.m_ion_mass)/E/grids[i_n];
    grids[divEcond] = grids[divE] - 4.*PI*grids[rho_c];

}

void Ideal2F::catchNullFieldDirection(std::vector<Grid> &grids) const
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

void Ideal2F::recomputeDT(std::vector<Grid>& grids) const
{
    // diagonal grid size
    const Grid& dx = m_pd.m_grids[PlasmaDomain::d_x];
    const Grid& dy = m_pd.m_grids[PlasmaDomain::d_y];
    Grid dr = (dx.square() + dy.square()).sqrt();
    // smallest wavenumber supported by grid structure
    Grid lambda = dr*2;
    Grid k = 2.*PI/lambda;
    // determine largest fluid velocity
    Grid v_mag_i = (grids[i_v_x].square() + grids[i_v_y].square()).sqrt();
    Grid v_mag_e = (grids[e_v_x].square() + grids[e_v_y].square()).sqrt();
    Grid v = v_mag_e.max(v_mag_i);
    // compute group velocity of Langmuir wave
    Grid v_th = (K_B*grids[e_temp]/M_ELECTRON).sqrt();
    Grid w_pe = (4*PI*grids[e_n]*E*E/M_ELECTRON).sqrt();
    Grid w_L = (3*k.square()*v_th.square()).sqrt();
    Grid v_L = w_L/k;
    // get timesteps
    Grid dt_wave = 1./w_pe;
    Grid dt_v = dr/(v + v_L);
    Grid dt_EM = dx*dy/(dx+dy)/C;
    grids[dt] = dt_wave.min(dt_v);
    if (!m_use_sub_cycling && !m_remove_curl_terms) grids[dt] = grids[dt].min(dt_EM);
    // record ion timescale
    grids[dt_i] = ionTimescale();
}

// returns timescale of bulk flow and acoustic waves for ions
Grid Ideal2F::ionTimescale() const
{
    Grid c_s = (m_pd.m_adiabatic_index*m_grids[i_press]/m_grids[i_rho]).sqrt();
    Grid vel_mag = (m_grids[i_v_x].square() + m_grids[i_v_y].square()).sqrt();
    Grid diagonals = (m_pd.m_grids[PlasmaDomain::d_x].square() + m_pd.m_grids[PlasmaDomain::d_y].square()).sqrt();
    return diagonals/(vel_mag + c_s);
}

std::vector<Grid> Ideal2F::subcycleMaxwell(const std::vector<Grid>& grids, const std::vector<Grid>& dj_tot, double step) const
{
    // determine sub-cycle timestep - Courant condition in two dimensions
    const Grid& dx = m_pd.m_grids[PlasmaDomain::d_x];
    const Grid& dy = m_pd.m_grids[PlasmaDomain::d_y];
    double dt_EM = m_epsilon_courant*(dx*dy/(dx+dy)/C).min(); // ideal timestep to satisfy Courant condition
    int num_steps = step/dt_EM + 1; // number of sub-cycles that satisfies Courant condition
    double dt = step/num_steps;
    #if VERBOSE
    if (m_verbose) std::cout << "Number Subcycles: " << num_steps << std::endl;
    #endif
    // preallocate variables
    std::vector<Grid> dEM_dt(6,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    std::vector<Grid> EM {grids[E_x],grids[E_y],grids[E_z],grids[b_x],grids[b_y],grids[b_z]};
    std::vector<Grid> EM_step(6,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    std::vector<Grid> EM_half = EM;
    std::vector<Grid> EM_laplacian(6,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    std::vector<Grid> j {grids[j_x],grids[j_y]};
    std::vector<Grid> dj {dj_tot[0]/num_steps,dj_tot[1]/num_steps};
    populate_boundary(dj[0]);
    populate_boundary(dj[1]);
    // loop over sub-cycles
    for (int i=0; i<num_steps; i++){
        // midpoint RK2 step
        maxwellCurlEqs(EM,j,EM_laplacian,dEM_dt);
        for (int k=0; k<EM.size(); k++) EM_half[k] = EM[k] + dEM_dt[k]*(dt/2.);
        maxwellCurlEqs(EM_half,j,EM_laplacian,dEM_dt);
        for (int k=0; k<EM.size(); k++){
            EM[k] += dEM_dt[k]*dt;
            EM_step[k] += dEM_dt[k]*dt;
        }
        // increment current density
        j[0] += dj[0];
        j[1] += dj[1];
    }
    // compute change in EM fields
    return EM_step;
}

void Ideal2F::maxwellCurlEqs(const std::vector<Grid>& EM,const std::vector<Grid>& j, std::vector<Grid>& EM_laplacian, std::vector<Grid>& dEM_dt) const
{
    dEM_dt[0] = + C*m_pd.derivative1D(EM[5],1) - 4.*PI*j[0];
    dEM_dt[1] = - C*m_pd.derivative1D(EM[5],0) - 4.*PI*j[1];
    dEM_dt[2] = + C*(m_pd.derivative1D(EM[4],0)-m_pd.derivative1D(EM[3],1));
    dEM_dt[3] = - C*m_pd.derivative1D(EM[2],1);
    dEM_dt[4] = + C*m_pd.derivative1D(EM[2],0);
    dEM_dt[5] = + C*(m_pd.derivative1D(EM[0],1) - m_pd.derivative1D(EM[1],0));
}

void Ideal2F::populate_boundary(Grid& grid) const
{
    // treat left boundary
    for (int i=0; i<m_pd.m_xl; i++){
        for (int j=0; j<grid.cols(); j++){
            grid(i,j) = grid(m_pd.m_xl,j);
        }
    }
    // treat right boundary
    for (int i=m_pd.m_xu+1; i<grid.rows(); i++){
        for (int j=0; j<grid.cols(); j++){
            grid(i,j) = grid(m_pd.m_xu,j);
        }
    }
    // treat bottom boundary
    for (int i=0; i<grid.rows(); i++){
        for (int j=0; j<m_pd.m_yl; j++){
            grid(i,j) = grid(i,m_pd.m_yl);
        }
    }
    // treat top boundary
    for (int i=0; i<grid.rows(); i++){
        for (int j=m_pd.m_yu+1; j<grid.cols(); j++){
            grid(i,j) = grid(i,m_pd.m_yu);
        }
    }
}
