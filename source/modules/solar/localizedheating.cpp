#include "localizedheating.hpp"
#include "plasmadomain.hpp"
#include "solarutils.hpp"
#include "idealmhd.hpp"
#include <sstream>
#include <iostream>

LocalizedHeating::LocalizedHeating(PlasmaDomain &pd) : Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    heating_template = Grid::Zero(1,1);
}

void LocalizedHeating::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "start_time") start_time = std::stod(this_rhs);
        else if(this_lhs == "duration") duration = std::stod(this_rhs);
        else if(this_lhs == "max_heating_rate") max_heating_rate = std::stod(this_rhs);
        else if(this_lhs == "stddev_x") stddev_x = std::stod(this_rhs);
        else if(this_lhs == "stddev_y") stddev_y = std::stod(this_rhs);
        else if(this_lhs == "center_x") center_x = std::stod(this_rhs);
        else if(this_lhs == "center_y") center_y = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void LocalizedHeating::preIterateModule(double dt){
    if(heating_template.size() != 1) return;
    heating_template = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*max_heating_rate,
                                                max_heating_rate,stddev_x,stddev_y,center_x,center_y).max(0.0);
    if(m_pd.x_bound_1 == PlasmaDomain::BoundaryCondition::Periodic){
        Grid left_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*max_heating_rate,
                                                max_heating_rate,stddev_x,stddev_y,center_x+(double)(m_pd.m_xdim),center_y).max(0.0);
        Grid right_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*max_heating_rate,
                                                max_heating_rate,stddev_x,stddev_y,center_x-(double)(m_pd.m_xdim),center_y).max(0.0);
        heating_template = heating_template.max(left_shifted).max(right_shifted);
    }
    if(m_pd.y_bound_1 == PlasmaDomain::BoundaryCondition::Periodic){
        Grid down_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*max_heating_rate,
                                                max_heating_rate,stddev_x,stddev_y,center_x,center_y+(double)(m_pd.m_ydim)).max(0.0);
        Grid up_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*max_heating_rate,
                                                max_heating_rate,stddev_x,stddev_y,center_x,center_y-(double)(m_pd.m_ydim)).max(0.0);
        heating_template = heating_template.max(down_shifted).max(up_shifted);
    }
}

void LocalizedHeating::postIterateModule(double dt){
    if(m_pd.m_time < start_time || m_pd.m_time > start_time+duration) return;
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) += m_pd.m_ghost_zone_mask*(dt*heating_template);
    m_pd.m_eqs->propagateChanges();
}

std::string LocalizedHeating::commandLineMessage() const
{
    std::ostringstream oss;
    oss.precision(4);
    oss << center_x << "," << center_y;
    std::string result = "Heating at " + oss.str();
    if(m_pd.m_time < start_time || m_pd.m_time > start_time+duration) result += " Off";
    else result += " On";
    return result;
}