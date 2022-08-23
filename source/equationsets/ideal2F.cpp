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
        else if (lhs[i] == "smooth_vars") m_smooth_vars = (rhs[i] == "true");
        else if (lhs[i] == "verbose_2F") m_verbose = (rhs[i] == "true");
        else if (lhs[i] == "apply_fixed_curl_bc") m_apply_fixed_curl_bc = (rhs[i] == "true");
        else{
            std::cerr << lhs[i] << " is not recognized for this equation set." << std::endl;
            assert(false);
        }
    }
}

void Ideal2F::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // evolve fluid variables
    grids[i_rho] += step*time_derivatives[0];
    grids[e_rho] += step*time_derivatives[1];
    grids[i_mom_x] += step*time_derivatives[2];
    grids[i_mom_y] += step*time_derivatives[3];
    grids[e_mom_x] += step*time_derivatives[4];
    grids[e_mom_y] += step*time_derivatives[5];
    grids[i_thermal_energy] += step*time_derivatives[6];
    grids[e_thermal_energy] += step*time_derivatives[7];
    // evolve electromagnetic fields
    if (m_use_sub_cycling){
        // compute change in current density
        Grid j_x_new = E*(grids[i_mom_x]/m_pd.m_ion_mass-grids[e_mom_x]/M_ELECTRON);
        Grid j_y_new = E*(grids[i_mom_y]/m_pd.m_ion_mass-grids[e_mom_y]/M_ELECTRON);
        std::vector<Grid> dj = {j_x_new - grids[j_x],j_y_new - grids[j_y]};
        // use sub-cycling to evolve E and B
        std::vector<Grid> step_EM = subcycleMaxwell(grids, dj, step);
        grids[E_x] += step_EM[0];
        grids[E_y] += step_EM[1];
        grids[E_z] += step_EM[2];
        grids[bi_x] += step_EM[3];
        grids[bi_y] += step_EM[4];
        grids[bi_z] += step_EM[5];
    }
    else{
        grids[E_x] += step*time_derivatives[8];
        grids[E_y] += step*time_derivatives[9];
        grids[E_z] += step*time_derivatives[10];
        grids[bi_x] += step*time_derivatives[11];
        grids[bi_y] += step*time_derivatives[12];
        grids[bi_z] += step*time_derivatives[13];
    }
    // propagate changes
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
    Grid i_F_x =  E*grids[i_n]*(grids[E_x] + i_v_cross_B[0]/C);
    Grid i_F_y =  E*grids[i_n]*(grids[E_y] + i_v_cross_B[1]/C);
    Grid e_F_x = -E*grids[e_n]*(grids[E_x] + e_v_cross_B[0]/C);
    Grid e_F_y = -E*grids[e_n]*(grids[E_y] + e_v_cross_B[1]/C); 
    // momentum equations   
    Grid i_d_mom_x_dt = - m_pd.transportDivergence2D(grids[i_mom_x], v_i)
                        - m_pd.derivative1D(grids[i_press], 0)
                        + grids[i_rho]*grids[grav_x] + i_viscous_force_x + i_F_x;
    Grid i_d_mom_y_dt = - m_pd.transportDivergence2D(grids[i_mom_y], v_i)
                        - m_pd.derivative1D(grids[i_press], 1)
                        + grids[i_rho]*grids[grav_y] + i_viscous_force_y + i_F_y;
    Grid e_d_mom_x_dt = - m_pd.transportDivergence2D(grids[e_mom_x], v_e)
                        - m_pd.derivative1D(grids[e_press], 0)
                        + grids[e_rho]*grids[grav_x] + e_viscous_force_x + e_F_x;
    Grid e_d_mom_y_dt = - m_pd.transportDivergence2D(grids[e_mom_y], v_e)
                        - m_pd.derivative1D(grids[e_press], 1)
                        + grids[e_rho]*grids[grav_y] + e_viscous_force_y + e_F_y;
    // energy equations
    Grid i_d_thermal_dt =   - m_pd.transportDivergence2D(grids[i_thermal_energy],v_i)
                            - grids[i_press]*m_pd.divergence2D(v_i);
    Grid e_d_thermal_dt =   - m_pd.transportDivergence2D(grids[e_thermal_energy],v_e)
                            - grids[e_press]*m_pd.divergence2D(v_e);
    // optionally evolve electromagnetic fields
    Grid dEx_dt, dEy_dt, dEz_dt, dBx_dt, dBy_dt, dBz_dt;
    if (!m_use_sub_cycling){
        dEx_dt = +C*m_pd.derivative1D(grids[b_z],1) - 4.*PI*grids[j_x];
        dEy_dt = -C*m_pd.derivative1D(grids[b_z],0) - 4.*PI*grids[j_y];
        dEz_dt = +C*(m_pd.derivative1D(grids[b_y],0) - m_pd.derivative1D(grids[b_x],1));
        dBx_dt = -C*m_pd.derivative1D(grids[E_z],1);
        dBy_dt = +C*m_pd.derivative1D(grids[E_z],0);
        dBz_dt = +C*(m_pd.derivative1D(grids[E_x],1) - m_pd.derivative1D(grids[E_y],0));
    }
    // return time derivatives
    return {i_d_rho_dt,e_d_rho_dt,i_d_mom_x_dt,i_d_mom_y_dt,e_d_mom_x_dt,e_d_mom_y_dt,i_d_thermal_dt,e_d_thermal_dt, dEx_dt, dEy_dt, dEz_dt, dBx_dt, dBy_dt, dBz_dt};
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
    if (m_smooth_vars) smooth_vars(grids);
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
    grids[dn] = grids[i_n] - grids[e_n];
    grids[divE] = m_pd.divergence2D({grids[E_x],grids[E_y]});
    grids[divB] = m_pd.divergence2D({grids[b_x],grids[b_y]});
    grids[curlE_z] = -(m_pd.derivative1D(grids[E_x],1) - m_pd.derivative1D(grids[E_y],0));
    if (m_apply_fixed_curl_bc) apply_fixed_curl_bc(grids[curlE_z]);
    grids[j_x_sg] = m_pd.m_sg.smooth_interior(grids[j_x]);
    grids[E_x_sg] = m_pd.m_sg.smooth_interior(grids[E_x]);
    grids[E_y_sg] = m_pd.m_sg.smooth_interior(grids[E_y]);
    populate_boundary(grids[j_x_sg]);
    populate_boundary(grids[E_x_sg]);
    populate_boundary(grids[E_y_sg]);
    grids[curlE_z_sg] = -(m_pd.derivative1D(grids[E_x_sg],1) - m_pd.derivative1D(grids[E_y_sg],0));
    if (m_apply_fixed_curl_bc) apply_fixed_curl_bc(grids[curlE_z_sg]);
    grids[B_z_sg] = m_pd.m_sg.smooth_interior(grids[b_z]);
    populate_boundary(grids[B_z_sg]);
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
    // diagonal grid size
    const Grid& dx = m_pd.m_grids[PlasmaDomain::d_x];
    const Grid& dy = m_pd.m_grids[PlasmaDomain::d_y];
    Grid dr = (dx.square() + dy.square()).sqrt();
    // smallest wavenumber supported by grid structure
    Grid lambda = dr*2;
    Grid k = 2.*PI/lambda;
    // determine largest fluid velocity
    Grid v_mag_i = (m_grids[i_v_x].square() + m_grids[i_v_y].square()).sqrt();
    Grid v_mag_e = (m_grids[e_v_x].square() + m_grids[e_v_y].square()).sqrt();
    Grid v = v_mag_e.max(v_mag_i);
    // compute group velocity of Langmuir wave
    Grid v_th = (K_B*m_grids[e_temp]/M_ELECTRON).sqrt();
    Grid w_pe = (4*PI*m_grids[e_n]*E*E/M_ELECTRON).sqrt();
    Grid w_L = (3*k.square()*v_th.square()).sqrt();
    Grid v_L = w_L/k;
    // get timesteps
    Grid dt_wave = 1./w_pe;
    Grid dt_v = dr/(v + v_L);
    Grid dt_EM = dx*dy/(dx+dy)/C;
    m_grids[dt] = dt_wave.min(dt_v);
    if (!m_use_sub_cycling) m_grids[dt] = m_grids[dt].min(dt_EM);
}

std::vector<Grid> Ideal2F::subcycleMaxwell(const std::vector<Grid>& grids, const std::vector<Grid>& dj_tot, double step)
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
    std::vector<Grid> j {grids[j_x],grids[j_y]};
    std::vector<Grid> dj {dj_tot[0]/num_steps,dj_tot[1]/num_steps};
    populate_boundary(dj[0]);
    populate_boundary(dj[1]);
    // loop over sub-cycles
    for (int i=0; i<num_steps; i++){
        // midpoint RK2 step
        maxwellCurlEqs(EM,j,dEM_dt);
        for (int k=0; k<EM.size(); k++) EM_half[k] = EM[k] + dEM_dt[k]*(dt/2.);
        maxwellCurlEqs(EM_half,j,dEM_dt);
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

void Ideal2F::maxwellCurlEqs(const std::vector<Grid>& EM,const std::vector<Grid>& j, std::vector<Grid>& dEM_dt)
{
    dEM_dt[0] =  C*m_pd.derivative1D(EM[5],1) - 4.*PI*j[0];
    dEM_dt[1] = -C*m_pd.derivative1D(EM[5],0) - 4.*PI*j[1];
    dEM_dt[2] =  C*(m_pd.derivative1D(EM[4],0)-m_pd.derivative1D(EM[3],1));
    if (m_apply_fixed_curl_bc) apply_fixed_curl_bc(dEM_dt[2]);
    dEM_dt[3] = -C*m_pd.derivative1D(EM[2],1);
    dEM_dt[4] =  C*m_pd.derivative1D(EM[2],0);
    dEM_dt[5] =  C*(m_pd.derivative1D(EM[0],1) - m_pd.derivative1D(EM[1],0));
    if (m_apply_fixed_curl_bc) apply_fixed_curl_bc(dEM_dt[5]);
}

void Ideal2F::apply_fixed_curl_bc(Grid& grid) const
{
    // determine number of interior cells to include in fixed boundary
    int num_int_x = 1;
    int num_int_y = 1;
    if (num_int_x < m_pd.m_sg.Bx()) num_int_x = m_pd.m_sg.Bx();
    if (num_int_y < m_pd.m_sg.By()) num_int_y = m_pd.m_sg.By();
    // treat left boundary
    for (int i=0; i<m_pd.m_xl+num_int_x; i++){
        #pragma omp parallel for
        for (int j=0; j<grid.cols(); j++){
            grid(i,j) = grid(m_pd.m_xl+num_int_x,j);
        }
    }
    // treat right boundary
    for (int i=m_pd.m_xu+1-num_int_x; i<grid.rows(); i++){
        #pragma omp parallel for
        for (int j=0; j<grid.cols(); j++){
            grid(i,j) = grid(m_pd.m_xu-num_int_x,j);
        }
    }
    // treat bottom boundary
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        #pragma omp parallel for
        for (int j=0; j<m_pd.m_yl+num_int_y; j++){
            grid(i,j) = grid(i,m_pd.m_yl+num_int_y);
        }
    }
    // treat top boundary
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=m_pd.m_yu+1-num_int_y; j<grid.cols(); j++){
            grid(i,j) = grid(i,m_pd.m_yu-num_int_y);
        }
    }
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

void Ideal2F::smooth_vars(std::vector<Grid> &grids) const
{
    int sg_ghost_xl = m_pd.m_xl + m_pd.m_sg.Bx();
    int sg_ghost_yl = m_pd.m_yl + m_pd.m_sg.By();
    int sg_ghost_xu = m_pd.m_xu - m_pd.m_sg.Bx();
    int sg_ghost_yu = m_pd.m_yu - m_pd.m_sg.By();
    for (auto ind : vars_to_smooth) grids[ind] = m_pd.m_sg.smooth(grids[ind],sg_ghost_xl,sg_ghost_xu,sg_ghost_yl,sg_ghost_yu);
}