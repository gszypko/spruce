#include "module.hpp"
#include "plasmadomain.hpp"
#include "fieldheating.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <iostream>

FieldHeating::FieldHeating(PlasmaDomain &pd): Module(pd) {}

void FieldHeating::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "rate") rate = std::stod(this_rhs);
        else if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "current_mode") current_mode = (this_rhs == "true");
        else std::cerr << this_lhs << " config not recognized.\n";
    }
    heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
}

void FieldHeating::computeHeating(){
    if(current_mode) heating = rate*C/(4.0*PI)*(m_pd.curl2D(m_pd.m_be_x+m_pd.grid(IdealMHD::bi_x),
                                                            m_pd.m_be_y+m_pd.grid(IdealMHD::bi_y))).abs();
    else heating = rate*(m_pd.grid(IdealMHD::b_magnitude)).square()/(8.0*PI);
}

void FieldHeating::postIterateModule(double dt){
    computeHeating();
    m_pd.grid(IdealMHD::thermal_energy) += m_pd.m_ghost_zone_mask*(dt*heating);
    m_pd.m_eqs->propagateChanges();
}

std::string FieldHeating::commandLineMessage() const
{
    return (current_mode ? "Field Heating On (Current Mode)" : "Field Heating On");
}

void FieldHeating::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const
{
    if(output_to_file){
        var_names.push_back("field_heating");
        var_grids.push_back(heating);
    }
}
