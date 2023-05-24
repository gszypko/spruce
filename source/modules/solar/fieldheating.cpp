#include "module.hpp"
#include "plasmadomain.hpp"
#include "fieldheating.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <iostream>
#include <cassert>

FieldHeating::FieldHeating(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    heating = Grid::Zero(1,1);
}

void FieldHeating::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "coeff") coeff = std::stod(this_rhs);
        else if(this_lhs == "current_pow") current_pow = std::stod(this_rhs);
        else if(this_lhs == "b_pow") b_pow = std::stod(this_rhs);
        else if(this_lhs == "n_pow") n_pow = std::stod(this_rhs);
        else if(this_lhs == "roc_pow") roc_pow = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void FieldHeating::setupModule(){
    heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
}

void FieldHeating::computeHeating(){
    if(coeff == 0.0) heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    else {
        heating = Grid(m_pd.m_xdim,m_pd.m_ydim,coeff);
        if(current_pow != 0.0) heating *= (C/(4.0*PI)*(m_pd.curl2D(m_pd.m_grids[PlasmaDomain::be_x]+m_pd.m_eqs->grid(IdealMHD::bi_x),
                                                                m_pd.m_grids[PlasmaDomain::be_y]+m_pd.m_eqs->grid(IdealMHD::bi_y))).abs()).pow(current_pow);
        if(b_pow != 0.0) heating *= (m_pd.m_eqs->grid(IdealMHD::b_mag)).pow(b_pow);
        if(n_pow != 0.0) heating *= (m_pd.m_eqs->grid(IdealMHD::n)).pow(n_pow);
        if(roc_pow != 0.0){
            Grid &b_hat_x = m_pd.m_eqs->grid(IdealMHD::b_hat_x), &b_hat_y = m_pd.m_eqs->grid(IdealMHD::b_hat_y);
            Grid curv_x = b_hat_x*m_pd.derivative1D(b_hat_x,0) + b_hat_y*m_pd.derivative1D(b_hat_x,1);
            Grid curv_y = b_hat_x*m_pd.derivative1D(b_hat_y,0) + b_hat_y*m_pd.derivative1D(b_hat_y,1);
            Grid roc = 1.0/((curv_x.square() + curv_y.square()).sqrt().max(1.0e-16));
            heating /= (roc).pow(roc_pow);
        }
    }
}

void FieldHeating::preIterateModule(double dt){
    computeHeating();
}

void FieldHeating::iterateModule(double dt){
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) += m_pd.m_ghost_zone_mask*(dt*heating);
    m_pd.m_eqs->propagateChanges();
    heating = m_pd.m_ghost_zone_mask*(dt*heating);
}


std::string FieldHeating::commandLineMessage() const
{
    return (coeff == 0.0 ? "Field Heating Zero" : "Field Heating On");
}

void FieldHeating::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if(output_to_file){
        if(heating.size() == 1) heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
        var_names.push_back("field_heating");
        var_grids.push_back(heating);
    }
}
