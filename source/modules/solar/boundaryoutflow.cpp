#include "boundaryoutflow.hpp"
#include "plasmadomain.hpp"
#include "solarutils.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <sstream>
#include <iostream>
#include <cmath>

BoundaryOutflow::BoundaryOutflow(PlasmaDomain &pd) : Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    accel_template = Grid::Zero(1,1);
}

void BoundaryOutflow::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "max_accel") max_accel = std::stod(this_rhs);
        else if(this_lhs == "falloff_length") falloff_length = std::stod(this_rhs);
        else if(this_lhs == "boundary") boundary = this_rhs;
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void BoundaryOutflow::setupModule(){
    accel_template = constructBoundaryAccel(max_accel,falloff_length);
}

void BoundaryOutflow::postIterateModule(double dt){
    if (boundary=="x_bound_1" || boundary=="x_bound_2")
        m_pd.m_eqs->grid(IdealMHD::mom_x) += m_pd.m_ghost_zone_mask*(dt*accel_template*m_pd.m_eqs->grid(IdealMHD::rho));
    else if (boundary=="y_bound_1" || boundary=="y_bound_2")
        m_pd.m_eqs->grid(IdealMHD::mom_y) += m_pd.m_ghost_zone_mask*(dt*accel_template*m_pd.m_eqs->grid(IdealMHD::rho));
    m_pd.m_eqs->propagateChanges();
}

std::string BoundaryOutflow::commandLineMessage() const
{
    return boundary + " boundary outflow enforced";
}

Grid BoundaryOutflow::constructBoundaryAccel(double strength,double length) const
{
    // initialize grids and references to PlasmaDomain grids
    const Grid& x = m_pd.grid("pos_x");
    const Grid& y = m_pd.grid("pos_y");

    if (boundary=="x_bound_1") return (-2.3*(x - x.min())/length).exp()*strength;
    else if (boundary=="x_bound_2") return (-2.3*(x - x.max())/length).exp()*strength;
    else if (boundary=="y_bound_2") return (-2.3*(y - y.max())/length).exp()*strength;
    else if (boundary=="y_bound_1") return (-2.3*(y - y.min())/length).exp()*strength;
    else assert(false && "BoundaryOutflow boundary config must be {x,y}_bound_{1,2}");
}
