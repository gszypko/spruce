#include "divcleaning.hpp"
#include "plasmadomain.hpp"
#include "idealmhd.hpp"
#include <iostream>
#include <cmath>

DivCleaning::DivCleaning(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
}

void DivCleaning::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "epsilon") epsilon = std::stod(this_rhs);
        if(this_lhs == "time_scale") time_scale = std::stod(this_rhs);
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void DivCleaning::setupModule(){
    const Grid &dx = m_pd.m_grids[PlasmaDomain::d_x], &dy = m_pd.m_grids[PlasmaDomain::d_y];
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);
    coeff = (one/(one/dx.square()+one/dy.square()))/2./time_scale;
}

void DivCleaning::postIterateModule(double dt){
    Grid b_x = m_pd.m_eqs->grid(IdealMHD::b_x), b_y = m_pd.m_eqs->grid(IdealMHD::b_y);
    Grid be_x = m_pd.m_grids[PlasmaDomain::be_x], be_y = m_pd.m_grids[PlasmaDomain::be_y];
    Grid bi_x = m_pd.m_eqs->grid(IdealMHD::bi_x), bi_y = m_pd.m_eqs->grid(IdealMHD::bi_y);
    int num_subcycles = computeNumSubcycles(dt);
    double dt_subcycle = dt/((double)num_subcycles);
    for(int subcycle = 0; subcycle < num_subcycles; subcycle++){
        bi_x += m_pd.m_ghost_zone_mask*dt_subcycle*coeff*(m_pd.derivative1D(m_pd.derivative1D(b_y,1),0)+m_pd.secondDerivative1D(b_x,0));
        bi_y += m_pd.m_ghost_zone_mask*dt_subcycle*coeff*(m_pd.derivative1D(m_pd.derivative1D(b_x,0),1)+m_pd.secondDerivative1D(b_y,1));
        b_x = bi_x + be_x;
        b_y = bi_y + be_y;
    }

    m_pd.m_eqs->grid(IdealMHD::bi_x) = bi_x;
    m_pd.m_eqs->grid(IdealMHD::bi_y) = bi_y;
    m_pd.m_eqs->propagateChanges();
}

int DivCleaning::computeNumSubcycles(double dt){
    return (int)(dt/(epsilon*time_scale)) + 1;
}

std::string DivCleaning::commandLineMessage() const
{
    return "Div. Cleaning On";
}
