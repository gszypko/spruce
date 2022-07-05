#include "eic_thermalization.hpp"

#define LENGTH_DECAY [](double r,double l){return exp(-r/l);} 

EICThermalization::EICThermalization(PlasmaDomain &pd): Module(pd) {}

void EICThermalization::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "eic_output_to_file") m_output_to_file = (rhs[i] == "true");
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
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

void EICThermalization::postIterateModule(double dt)
{
    const Grid& n = m_pd.m_eqs->grid("n");
    const Grid& Te = m_pd.m_eqs->grid("temp_e");
    Grid& eps_e = m_pd.m_eqs->grid("thermal_energy_e");
    Grid& eps_i = m_pd.m_eqs->grid("thermal_energy_i");

    m_vars[a] = ((3./4./PI)/n).pow(1./3.);
    m_vars[w_pe] = ((4.*PI*E*E/M_ELECTRON)*n).sqrt();
    m_vars[Gam_e] = (E*E/K_B)/Te/m_vars[a];
    m_vars[Lam_e] = (1./sqrt(3.))/m_vars[Gam_e].pow(3./2.);
    m_vars[gam_ei] = sqrt(2./3./PI)*m_vars[Gam_e].pow(3./2.)*m_vars[w_pe]*m_vars[Lam_e].log();
    m_vars[nu_ei] = (2.*M_ELECTRON/m_pd.m_ion_mass)*m_vars[gam_ei];
    m_vars[dEdt] = m_vars[nu_ei]*(eps_e-eps_i);
    m_vars[dEdEi] = m_vars[dEdt]*dt/eps_i;
    eps_e += -m_vars[dEdt]*dt;
    eps_i += m_vars[dEdt]*dt;
    m_pd.m_eqs->propagateChanges();
}

void EICThermalization::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if (m_output_to_file) {
        var_names.push_back(m_var_names[dEdEi]);
        var_grids.push_back(m_vars[dEdEi]);
    }
}