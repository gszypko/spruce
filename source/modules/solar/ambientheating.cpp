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
        else if(this_lhs == "split_exp_mode") split_exp_mode = (this_rhs == "true");
        else if(this_lhs == "split_exp_scale_height") split_exp_scale_height = std::stod(this_rhs);
        else if(this_lhs == "split_exp_start_height") split_exp_start_height = std::stod(this_rhs);
        else if(this_lhs == "ms_electron_heating_fraction") ms_electron_heating_fraction = std::stod(this_rhs);
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void AmbientHeating::setupModule(){
    if(exp_mode){
        heating = m_pd.m_ghost_zone_mask*exp_base_heating_rate*(-1.0*m_pd.m_grids[PlasmaDomain::pos_y]/exp_scale_height).exp();
        if(split_exp_mode){
            assert(split_exp_scale_height >= exp_scale_height && "Split (upper) scale height must be larger than lower one.");
            double split_base_heating_rate = exp_base_heating_rate*std::exp((exp_scale_height - split_exp_scale_height)*split_exp_start_height/(exp_scale_height*split_exp_scale_height));
            Grid other_heating = m_pd.m_ghost_zone_mask*split_base_heating_rate*(-1.0*m_pd.m_grids[PlasmaDomain::pos_y]/split_exp_scale_height).exp();
            heating = heating.max(other_heating);
        }
    }
    else heating = m_pd.m_ghost_zone_mask*heating_rate;
    assert(ms_electron_heating_fraction >= 0.0 && ms_electron_heating_fraction <= 1.0 && "Ambient Heating MS electron heating fraction must be between 0 and 1");
}

void AmbientHeating::postIterateModule(double dt){
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) += dt*heating;
    m_pd.m_eqs->propagateChanges();
    if(m_pd.m_multispecies_mode) {
        if(ms_electron_heating_fraction < 1.0) m_pd.m_cumulative_ion_heating += (1.0 - ms_electron_heating_fraction)*m_pd.m_ghost_zone_mask*heating*dt;
        if(ms_electron_heating_fraction > 0.0) m_pd.m_cumulative_electron_heating += ms_electron_heating_fraction*m_pd.m_ghost_zone_mask*heating*dt;
    }
}

std::string AmbientHeating::commandLineMessage() const
{
    return exp_mode ? "Ambient Heating On (Exp. Mode)" : "Ambient Heating On";
}
