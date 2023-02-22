#include "eic_thermalization.hpp"

#define LENGTH_DECAY [](double r,double l){return exp(-r/l);} 

EICThermalization::EICThermalization(PlasmaDomain &pd): Module(pd) {}

void EICThermalization::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    return;
}

void EICThermalization::setupModule()
{
    // initialize grids for module
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

void EICThermalization::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt)
{
    const Grid& n = grids[m_pd.m_eqs->name2index("n")];
    const Grid& Te = grids[m_pd.m_eqs->name2index("e_temp")];
    const Grid& eps_e = grids[m_pd.m_eqs->name2index("e_thermal_energy")];
    const Grid& eps_i = grids[m_pd.m_eqs->name2index("i_thermal_energy")];

    m_vars[a] = ((3./4./PI)/n).pow(1./3.);
    m_vars[w_pe] = ((4.*PI*E*E/M_ELECTRON)*n).sqrt();
    m_vars[Gam_e] = (E*E/K_B)/Te/m_vars[a];
    m_vars[Lam_e] = (1./sqrt(3.))/m_vars[Gam_e].pow(3./2.);
    m_vars[gam_ei] = sqrt(2./3./PI)*m_vars[Gam_e].pow(3./2.)*m_vars[w_pe]*m_vars[Lam_e].log();
    m_vars[nu_ei] = (2.*M_ELECTRON/m_pd.m_ion_mass)*m_vars[gam_ei];
    m_vars[dEdt] = m_vars[nu_ei]*(eps_e-eps_i);
    
    grids_dt[m_pd.m_eqs->name2evolvedindex("e_thermal_energy")] -= m_vars[dEdt]*m_pd.m_ghost_zone_mask;
    grids_dt[m_pd.m_eqs->name2evolvedindex("i_thermal_energy")] += m_vars[dEdt]*m_pd.m_ghost_zone_mask;
}