#include "coulomb_explosion.hpp"

#define LENGTH_DECAY [](double r,double l){return exp(-r/l);} 

CoulombExplosion::CoulombExplosion(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
}

void CoulombExplosion::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "timescale") m_timescale = std::stod(rhs[i]);
        else if(lhs[i] == "lengthscale") m_lengthscale = std::stod(rhs[i]);
        else if(lhs[i] == "strength") m_strength = std::stod(rhs[i]);
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void CoulombExplosion::setupModule()
{
    m_pd.grid(IdealMHD::Fcx) = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    m_pd.grid(IdealMHD::Fcy) = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    m_pd.grid(IdealMHD::dPx) = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    m_pd.grid(IdealMHD::dPy) = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
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
    return (4.*PI)*Grid::trapzcum(r,r.square()*rho_c);
}

void CoulombExplosion::postIterateModule(double dt)
{
    // references to grids
    const Grid& x {m_pd.m_internal_grids[PlasmaDomain::pos_x]};
    const Grid& y {m_pd.m_internal_grids[PlasmaDomain::pos_y]};
    const Grid& n {m_pd.grid(IdealMHD::n)};
    // compute distance from plasma center
    Grid r_sq = x.square()+y.square();
    Grid r = r_sq.sqrt();
    // bin the distances and densities
    Grid r_bin;
    Grid n_bin = Grid::bin_as_list(r,n,101,r_bin);
    // compute the charge density
    m_pd.grid(IdealMHD::dPx) = m_pd.derivative1D(m_pd.grid(IdealMHD::press), 0);
    m_pd.grid(IdealMHD::dPy) = m_pd.derivative1D(m_pd.grid(IdealMHD::press), 1);
    Grid rho_c = compute_charge_density(r_bin,n_bin);
    Grid Q_vec = compute_total_charge(r_bin,rho_c);
    // compute electric field
    Grid Q = Grid::interp_as_list(r_bin,Q_vec,r);
    Grid r_cubed = r_sq*r;
    Grid Ex = x*Q/r_cubed;
    Grid Ey = y*Q/r_cubed;
    // compute force profile
    m_pd.grid(IdealMHD::Fcx) = E*n*Ex;
    m_pd.grid(IdealMHD::Fcy) = E*n*Ey;
    // add forces
    m_pd.grid(IdealMHD::mom_x) += m_pd.grid(IdealMHD::Fcx)*dt;
    m_pd.grid(IdealMHD::mom_y) += m_pd.grid(IdealMHD::Fcy)*dt;
    m_pd.m_eqs->propagateChanges();
}

std::string CoulombExplosion::commandLineMessage() const
{
    return "Coulomb Explosion On";
}