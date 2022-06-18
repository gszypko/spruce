#include "coulomb_explosion.hpp"

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

// returns rho_c values in cgs at each distance in <r>
Grid CoulombExplosion::compute_charge_density(const Grid& r)
{
    Grid test;
    return test;
}

void CoulombExplosion::postIterateModule(double dt)
{
    // compute charge density
    const Grid& x {m_pd.m_internal_grids[PlasmaDomain::pos_x]};
    const Grid& y {m_pd.m_internal_grids[PlasmaDomain::pos_y]};
    const Grid& n {m_pd.grid(IdealMHD::n)};
    Grid r_sq = x.square()+y.square();
    Grid r = r_sq.sqrt();
    Grid r_bin = Grid::Linspace(r.min(),r.max(),101);
    double dr = r_bin(1)-r_bin(0);
    Grid edges = Grid::Linspace(r_bin.min()-dr/2.,r_bin.max()+dr/2.,r_bin.size()+1);
    Grid n_bin = Grid::bin_as_list(r,n,edges);
    m_pd.grid(IdealMHD::dPx) = m_pd.derivative1D(m_pd.grid(IdealMHD::press), 0);
    m_pd.grid(IdealMHD::dPy) = m_pd.derivative1D(m_pd.grid(IdealMHD::press), 1);

    m_pd.m_eqs->propagateChanges();
}

std::string CoulombExplosion::commandLineMessage() const
{
    return "Coulomb Explosion On";
}