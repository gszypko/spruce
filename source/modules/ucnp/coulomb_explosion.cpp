#include "coulomb_explosion.hpp"

#define LENGTH_DECAY [](double r,double l){return exp(-r/l);} 

CoulombExplosion::CoulombExplosion(PlasmaDomain &pd): Module(pd) {}

void CoulombExplosion::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "timescale") m_timescale = std::stod(rhs[i]);
        else if(lhs[i] == "lengthscale") m_lengthscale = std::stod(rhs[i]);
        else if(lhs[i] == "strength") m_strength = std::stod(rhs[i]);
        else if(lhs[i] == "output_to_file") output_to_file = (rhs[i] == "true");
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void CoulombExplosion::setupModule()
{
    // initialize internal grids
    m_vars.resize(num_vars,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    // ensure that grid dependencies exist within EquationSet
    std::vector<std::string> grid_names = m_pd.m_eqs->allNames();
    for (auto name : m_eqset_grids){
        auto it = std::find(grid_names.begin(),grid_names.end(),name);
        if (it==grid_names.end()){
            std::cerr << "Grid <" << name << "> was not found within the EquationSet." << std::endl;
            assert(false);
        }
    }
}

// returns rho_c values in cgs at each distance in <r>
Grid CoulombExplosion::compute_charge_density(const Grid& r, const Grid& n) const
{
    double time_decay = exp(-m_pd.m_time/m_timescale);
    Grid length_decay = r.for_each(m_lengthscale,[](double a,double b){return exp(-a/b);});
    Grid rho_c = (m_strength*time_decay*E)*n*length_decay;
    return rho_c;
}

Grid CoulombExplosion::compute_total_charge(const Grid& r,const Grid& rho_c) const
{
    return (4.*PI)*Grid::TrapzCum(r,r.square()*rho_c);
}

void CoulombExplosion::postIterateModule(double dt)
{
    if (m_pd.m_time < 3*m_timescale){
    // references to 2D grids
    const Grid& x {m_pd.m_grids[PlasmaDomain::pos_x]};
    const Grid& y {m_pd.m_grids[PlasmaDomain::pos_y]};
    const Grid& n {m_pd.m_eqs->grid("n")};
    const Grid& press {m_pd.m_eqs->grid("press")};
    Grid& mom_x {m_pd.m_eqs->grid("mom_x")};
    Grid& mom_y {m_pd.m_eqs->grid("mom_y")};
    // compute distance from plasma center
    Grid r_sq = x.square()+y.square();
    Grid r = r_sq.sqrt();
    // bin the distances and densities
    Grid r_bin;
    Grid n_bin = Grid::Bin1D(r,n,101,r_bin);
    // compute the charge density
    m_vars[dP_x] = m_pd.derivative1D(press, 0);
    m_vars[dP_y] = m_pd.derivative1D(press, 1);
    Grid rho_c_vec = compute_charge_density(r_bin,n_bin);
    Grid Q_vec = compute_total_charge(r_bin,rho_c_vec);
    // compute electric field
    Grid Q = Grid::Interp1D(r_bin,Q_vec,r);
    Grid rho_c = Grid::Interp1D(r_bin,rho_c_vec,r);
    Grid r_cubed = r_sq*r;
    Grid Ex = x*Q/r_cubed;
    Grid Ey = y*Q/r_cubed;
    // compute force profile
    m_vars[F_x] = rho_c*Ex;
    m_vars[F_y] = rho_c*Ey;
    // add forces
    mom_x += m_vars[F_x]*dt;
    mom_y += m_vars[F_y]*dt;
    m_pd.m_eqs->propagateChanges();
    }
}

std::string CoulombExplosion::commandLineMessage() const
{
    return "Coulomb Explosion On";
}

void CoulombExplosion::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if (output_to_file) {
        for (int i=0; i<num_vars; i++){
            var_names.push_back(m_var_names[i]);
            var_grids.push_back(m_vars[i]);
        }
    }
}