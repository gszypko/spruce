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
        m_pd.m_eqs->grid(IdealMHD::mom_x) += (dt*accel_template*m_pd.m_eqs->grid(IdealMHD::rho));
    else if (boundary=="y_bound_1" || boundary=="y_bound_2")
        m_pd.m_eqs->grid(IdealMHD::mom_y) += (dt*accel_template*m_pd.m_eqs->grid(IdealMHD::rho));
    m_pd.m_eqs->propagateChanges();
}

std::string BoundaryOutflow::commandLineMessage() const
{
    return boundary + " boundary outflow enforced";
}

Grid BoundaryOutflow::constructBoundaryAccel(double strength,double length) const
{
    const Grid& x = m_pd.grid("pos_x");
    const Grid& y = m_pd.grid("pos_y");

    Grid ghost_mask_withboundary = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    int xl = m_pd.m_xl, xu = m_pd.m_xu, yl = m_pd.m_yl, yu = m_pd.m_yu;
    if (boundary=="x_bound_1") xl = m_pd.m_xl_dt;
    else if (boundary=="x_bound_2") xu = m_pd.m_xu_dt;
    else if (boundary=="y_bound_2") yu = m_pd.m_yu_dt;
    else if (boundary=="y_bound_1") yl = m_pd.m_yl_dt;
    for(int i = xl; i <= xu; i++){
        for(int j = yl; j <= yu; j++){
        ghost_mask_withboundary(i,j) = 1.0;
        }
    }

    Grid result;
    if (boundary=="x_bound_1") result = (-2.3*(x - x.min())/length).exp()*strength;
    else if (boundary=="x_bound_2") result = (2.3*(x - x.max())/length).exp()*strength;
    else if (boundary=="y_bound_2") result = (2.3*(y - y.max())/length).exp()*strength;
    else if (boundary=="y_bound_1") result = (-2.3*(y - y.min())/length).exp()*strength;
    else assert(false && "BoundaryOutflow boundary config must be {x,y}_bound_{1,2}");
    return ghost_mask_withboundary*result;
}
