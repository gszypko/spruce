#include "momentuminjection.hpp"
#include "plasmadomain.hpp"
#include "solarutils.hpp"
#include "idealmhd.hpp"
#include <sstream>
#include <iostream>
#include <cmath>
MomentumInjection::MomentumInjection(PlasmaDomain &pd) : Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    accel_template = std::vector<Grid>(2,Grid::Zero(1,1));
    dir = std::vector<double>(2,0.0);
}

void MomentumInjection::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "start_time") start_time = std::stod(this_rhs);
        else if(this_lhs == "duration") duration = std::stod(this_rhs);
        else if(this_lhs == "max_accel") max_accel = std::stod(this_rhs);
        else if(this_lhs == "stddev_x") stddev_x = std::stod(this_rhs);
        else if(this_lhs == "stddev_y") stddev_y = std::stod(this_rhs);
        else if(this_lhs == "center_x") center_x = std::stod(this_rhs);
        else if(this_lhs == "center_y") center_y = std::stod(this_rhs);
        else if(this_lhs == "dir_x") dir[0] = std::stod(this_rhs);
        else if(this_lhs == "dir_y") dir[1] = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void MomentumInjection::setupModule(){
    assert(!(dir[0] == 0.0 && dir[1] == 0.0) && "Momentum Injection module must be given a nonzero acceleration direction");
    //normalize direction vector
    double mag = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1]);
    dir[0] = dir[0] / mag;
    dir[1] = dir[1] / mag;
    //compute acceleration template
    for(int i : {0,1}){
        double curr_max_accel = dir[i]*max_accel;
        accel_template[i] = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*curr_max_accel,
                                                    curr_max_accel,stddev_x,stddev_y,center_x,center_y).max(0.0);
        if(m_pd.x_bound_1 == PlasmaDomain::BoundaryCondition::Periodic){
            Grid left_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*curr_max_accel,
                                                    curr_max_accel,stddev_x,stddev_y,center_x+(double)(m_pd.m_xdim),center_y).max(0.0);
            Grid right_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*curr_max_accel,
                                                    curr_max_accel,stddev_x,stddev_y,center_x-(double)(m_pd.m_xdim),center_y).max(0.0);
            accel_template[i] = accel_template[i].max(left_shifted).max(right_shifted);
        }
        if(m_pd.y_bound_1 == PlasmaDomain::BoundaryCondition::Periodic){
            Grid down_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*curr_max_accel,
                                                    curr_max_accel,stddev_x,stddev_y,center_x,center_y+(double)(m_pd.m_ydim)).max(0.0);
            Grid up_shifted = SolarUtils::GaussianGrid(m_pd.m_xdim,m_pd.m_ydim,-1.0e-2*curr_max_accel,
                                                    curr_max_accel,stddev_x,stddev_y,center_x,center_y-(double)(m_pd.m_ydim)).max(0.0);
            accel_template[i] = accel_template[i].max(down_shifted).max(up_shifted);
        }
    }
}

void MomentumInjection::postIterateModule(double dt){
    if(m_pd.m_time < start_time || m_pd.m_time > start_time+duration) return;
    m_pd.m_eqs->grid(IdealMHD::mom_x) += m_pd.m_ghost_zone_mask*(dt*accel_template[0]*m_pd.m_eqs->grid(IdealMHD::rho));
    m_pd.m_eqs->grid(IdealMHD::mom_y) += m_pd.m_ghost_zone_mask*(dt*accel_template[1]*m_pd.m_eqs->grid(IdealMHD::rho));
    m_pd.m_eqs->propagateChanges();
}

std::string MomentumInjection::commandLineMessage() const
{
    std::ostringstream oss;
    oss.precision(4);
    oss << center_x << "," << center_y;
    std::string result = "Momentum injection at " + oss.str();
    if(m_pd.m_time < start_time || m_pd.m_time > start_time+duration) result += " Off";
    else result += " On";
    return result;
}