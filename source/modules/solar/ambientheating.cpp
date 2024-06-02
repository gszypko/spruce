#include "ambientheating.hpp"
#include "plasmadomain.hpp"
#include "idealmhd.hpp"
#include <iostream>
#include <cmath>

AmbientHeating::AmbientHeating(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
}

void AmbientHeating::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "heating_rate") heating_rate = std::stod(this_rhs);
        if(this_lhs == "exp_mode") exp_mode = (this_rhs == "true");
        else if(this_lhs == "exp_base_heating_rate") exp_base_heating_rate = std::stod(this_rhs);
        else if(this_lhs == "exp_scale_height") exp_scale_height = std::stod(this_rhs);
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void AmbientHeating::setupModule(){
    if(exp_mode) heating = m_pd.m_ghost_zone_mask*exp_base_heating_rate*(-1.0*m_pd.m_grids[PlasmaDomain::pos_y]/exp_scale_height).exp();
    else heating = m_pd.m_ghost_zone_mask*heating_rate;
}

void AmbientHeating::postIterateModule(double dt){
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) += dt*heating;
    m_pd.m_eqs->propagateChanges();
}

std::string AmbientHeating::commandLineMessage() const
{
    return exp_mode ? "Ambient Heating On (Exp. Mode)" : "Ambient Heating On";
}
