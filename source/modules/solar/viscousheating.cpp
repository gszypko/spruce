#include "module.hpp"
#include "plasmadomain.hpp"
#include "viscousheating.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <iostream>
#include <cassert>

ViscousHeating::ViscousHeating(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    heating = Grid::Zero(1,1);
}

void ViscousHeating::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "coeff") coeff = std::stod(this_rhs);
        else if(this_lhs == "inactive_mode") inactive_mode = (this_rhs == "true");
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void ViscousHeating::setupModule(){
    heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
}

void ViscousHeating::computeHeating(){
    if(coeff == 0.0) heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    else {
        Grid &b_hat_x = m_pd.m_eqs->grid(IdealMHD::b_hat_x), &b_hat_y = m_pd.m_eqs->grid(IdealMHD::b_hat_y), &b_hat_z = m_pd.m_eqs->grid(IdealMHD::b_hat_z);        
        Grid &v_x = m_pd.m_eqs->grid(IdealMHD::v_x), &v_y = m_pd.m_eqs->grid(IdealMHD::v_y), &v_z = m_pd.m_eqs->grid(IdealMHD::v_z);
        Grid del_x_v_x = m_pd.derivative1D(v_x,0), del_x_v_y = m_pd.derivative1D(v_y,0), del_x_v_z = m_pd.derivative1D(v_z,0),
             del_y_v_x = m_pd.derivative1D(v_x,1), del_y_v_y = m_pd.derivative1D(v_y,1), del_y_v_z = m_pd.derivative1D(v_z,1);
        Grid scalar = -3.0*coeff*m_pd.m_eqs->grid(IdealMHD::temp).pow(2.5)*(
            b_hat_x*(b_hat_x*del_x_v_x + b_hat_y*del_y_v_x)
            + b_hat_y*(b_hat_x*del_x_v_y + b_hat_y*del_y_v_y)
            + b_hat_z*(b_hat_x*del_x_v_z + b_hat_y*del_y_v_z)
            - (del_x_v_x + del_y_v_y)/3.0
        );
        heating = scalar*(
            (1.0/3.0 - b_hat_x.square())*del_x_v_x + (1.0/3.0 - b_hat_y.square())*del_y_v_y
            -b_hat_x*b_hat_y*del_y_v_x - b_hat_y*b_hat_x*del_x_v_y - b_hat_z*b_hat_x*del_x_v_z - b_hat_z*b_hat_y*del_y_v_z
        );
    }
}

void ViscousHeating::preIterateModule(double dt){
    computeHeating();
}

void ViscousHeating::iterateModule(double dt){
    heating = m_pd.m_ghost_zone_mask*(dt*heating);
    if(inactive_mode) return;
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) += heating;
    m_pd.m_eqs->propagateChanges();
    if(m_pd.m_multispecies_mode) m_pd.m_cumulative_ion_heating += heating;
}


std::string ViscousHeating::commandLineMessage() const
{
    std::string message = "Viscous Heating";
    if(coeff == 0.0) message += " Zero";
    else message += " On";
    if(inactive_mode) message += " (Not Applied)";
    return message;
}

void ViscousHeating::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if(output_to_file){
        if(heating.size() == 1) heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
        var_names.push_back("viscous_heating");
        var_grids.push_back(heating);
    }
}
