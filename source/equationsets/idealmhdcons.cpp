#include "idealmhdcons.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "grid.hpp"
#include <vector>
#include <iostream>

IdealMHDCons::IdealMHDCons(PlasmaDomain &pd): EquationSet(pd,def_var_names(),evolved_var_names()) {}

void IdealMHDCons::parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for (int i=0; i<lhs.size(); i++){
        if (lhs[i] == "global_viscosity") m_global_viscosity = stod(rhs[i]);
        else{
            std::cerr << lhs[i] << " is not recognized for this equation set." << std::endl;
            assert(false);
        }
    }
}

void IdealMHDCons::computeTimeDerivativesDerived(const std::vector<Grid> &grids, std::vector<Grid> &grids_dt){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    // PlasmaDomain grid references for more concise notation
    Grid& d_x = m_pd.m_grids[PlasmaDomain::d_x];
    Grid& d_y = m_pd.m_grids[PlasmaDomain::d_y];
    Grid& be_x = m_pd.m_grids[PlasmaDomain::be_x];
    Grid& be_y = m_pd.m_grids[PlasmaDomain::be_y];
    // continuity equations
    std::vector<Grid> v = {grids[v_x],grids[v_y]};
    grids_dt[ev2i("rho")] = -m_pd.transportDivergence2D(grids[rho],v);
    // tensor for momentum equation
    Grid Txx = grids[press_tot] + grids[b_x]*grids[b_x]/(4.*PI);
    Grid Tyy = grids[press_tot] + grids[b_y]*grids[b_y]/(4.*PI);
    Grid Txy = grids[b_x]*grids[b_y]/(4.*PI);
    Grid Tyx = grids[b_y]*grids[b_x]/(4.*PI);
    // momentum equations   
    grids_dt[ev2i("mom_x")] =   - m_pd.transportDivergence2D(grids[mom_x], v)
                        + m_pd.m_ghost_zone_mask * (grids[rho]*grids[grav_x] - m_pd.divergence2D({Txx,Tyx}));
    grids_dt[ev2i("mom_y")] =   - m_pd.transportDivergence2D(grids[mom_y], v)
                        + m_pd.m_ghost_zone_mask * (grids[rho]*grids[grav_y] - m_pd.divergence2D({Txy,Tyy}));
    // energy equations
    Grid v_dot_B = Grid::DotProduct2D(v,{grids[b_x],grids[b_y]});
    Grid Ux = v_dot_B*grids[b_x]/(4.*PI)-grids[press_tot]*grids[v_x];
    Grid Uy = v_dot_B*grids[b_y]/(4.*PI)-grids[press_tot]*grids[v_y];
    grids_dt[ev2i("energy")] =  - m_pd.transportDivergence2D(grids[energy],v)
                        + m_pd.divergence2D({Ux,Uy});
    // induction equations
    Grid Yxx = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid Yyy = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    Grid Yxy = grids[v_x]*grids[b_y] - grids[b_x]*grids[v_y];
    Grid Yyx = -1.*Yxy;
    grids_dt[ev2i("bi_x")] = -m_pd.divergence2D({Yxx,Yyx});
    grids_dt[ev2i("bi_y")] = -m_pd.divergence2D({Yyy,Yxy});

    // apply ghost zone mask
    for (int i=0; i<evolved_var_names().size(); i++) grids_dt[i] *= m_pd.m_ghost_zone_mask;

}

void IdealMHDCons::propagateChanges(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    m_pd.updateGhostZones();
    recomputeDerivedVarsFromEvolvedVars(grids);
    recomputeDT();
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

void IdealMHDCons::catchNullFieldDirection(std::vector<Grid> &grids)
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

void IdealMHDCons::recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[press] = 2*grids[n]*K_B*grids[temp].max(m_pd.temp_min);
    grids[thermal_energy] = grids[press]/(m_pd.m_adiabatic_index - 1.0);
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[kinetic_energy] = 0.5*grids[rho]*(grids[v_x].square()+grids[v_y].square());
    grids[b_x] = m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x];
    grids[b_y] = m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_mag] = (grids[b_x].square() + grids[b_y].square()).sqrt();
    grids[mag_energy] = grids[b_mag].square()/(8.*PI);
    grids[energy] = grids[thermal_energy] + grids[kinetic_energy] + grids[mag_energy];
}

void IdealMHDCons::recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    grids[n] = (grids[rho]/m_pd.m_ion_mass).max(m_pd.density_min);
    grids[rho] = grids[n]*m_pd.m_ion_mass;
    grids[v_x] = grids[mom_x]/grids[rho];
    grids[v_y] = grids[mom_y]/grids[rho];
    grids[kinetic_energy] = 0.5*grids[rho]*(grids[v_x].square()+grids[v_y].square());
    grids[b_x] = m_pd.m_grids[PlasmaDomain::be_x] + grids[bi_x];
    grids[b_y] = m_pd.m_grids[PlasmaDomain::be_y] + grids[bi_y];
    grids[b_mag] = (grids[b_x].square() + grids[b_y].square()).sqrt();
    grids[mag_energy] = grids[b_mag].square()/(8.*PI);
    grids[b_hat_x] = grids[b_x]/grids[b_mag];
    grids[b_hat_y] = grids[b_y]/grids[b_mag];
    grids[thermal_energy] = grids[energy] - grids[mag_energy] - grids[kinetic_energy].max(m_pd.thermal_energy_min);
    grids[press] = (m_pd.m_adiabatic_index - 1.0)*grids[thermal_energy];
    grids[temp] = (grids[press]/(2*K_B*grids[n])).max(m_pd.temp_min);
    catchNullFieldDirection(grids);
}