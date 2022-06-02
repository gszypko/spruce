#include "ambientheating.hpp"
#include "plasmadomain.hpp"
#include "idealmhd.hpp"
#include <iostream>

AmbientHeating::AmbientHeating(PlasmaDomain &pd): Module(pd) {}

void AmbientHeating::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "heating_rate") heating_rate = std::stod(rhs[i]);
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void AmbientHeating::postIterateModule(double dt){
    m_pd.grid(IdealMHD::thermal_energy) += m_pd.m_ghost_zone_mask*(dt*heating_rate);
    m_pd.m_eqs->propagateChanges();
}

std::string AmbientHeating::commandLineMessage() const
{
    return "Ambient Heating On";
}