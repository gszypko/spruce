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
        else if(this_lhs == "rad_mirror_mode") rad_mirror_mode = (this_rhs == "true");
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void AmbientHeating::setupModule(){
    if(rad_mirror_mode) heating = m_pd.m_ghost_zone_mask*radiativeLosses();
    else heating = m_pd.m_ghost_zone_mask*heating_rate;
}

void AmbientHeating::postIterateModule(double dt){
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) += dt*heating;
    m_pd.m_eqs->propagateChanges();
}

std::string AmbientHeating::commandLineMessage() const
{
    return rad_mirror_mode ? "Ambient Heating On (Matching Initial Rad. Loss)" : "Ambient Heating On";
}

Grid AmbientHeating::radiativeLosses(){
    double cutoff_temp = 3.0e4;
    double cutoff_ramp = 1.0e3;
    Grid result = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    #pragma omp parallel for collapse(2)
    for (int i = m_pd.m_xl; i <= m_pd.m_xu; i++){
        for(int j = m_pd.m_yl; j <= m_pd.m_yu; j++){
            if(m_pd.m_eqs->grid(IdealMHD::temp)(i,j) < cutoff_temp) result(i,j) = 0.0;
            else {
                double logtemp = std::log10(m_pd.m_eqs->grid(IdealMHD::temp)(i,j));
                double n = m_pd.m_eqs->grid(IdealMHD::rho)(i,j)/m_pd.m_ion_mass;
                double chi, alpha;
                if(logtemp <= 4.97){
                    chi = 1.09e-31;
                    alpha = 2.0;
                } else if(logtemp <= 5.67){
                    chi = 8.87e-17;
                    alpha = -1.0;
                } else if(logtemp <= 6.18){
                    chi = 1.90e-22;
                    alpha = 0.0;
                } else if(logtemp <= 6.55){
                    chi = 3.53e-13;
                    alpha = -1.5;
                } else if(logtemp <= 6.90){
                    chi = 3.46e-25;
                    alpha = 1.0/3.0;
                } else if(logtemp <= 7.63){
                    chi = 5.49e-16;
                    alpha = -1.0;
                } else{
                    chi = 1.96e-27;
                    alpha = 0.5;
                }
                result(i,j) = n*n*chi*std::pow(m_pd.m_eqs->grid(IdealMHD::temp)(i,j),alpha);
                if(m_pd.m_eqs->grid(IdealMHD::temp)(i,j) < cutoff_temp + cutoff_ramp){
                    double ramp = (m_pd.m_eqs->grid(IdealMHD::temp)(i,j) - cutoff_temp)/cutoff_ramp;
                    result(i,j) *= ramp;
                }
            }
        }
    }
    return result;
}
