#include "viscosity.hpp"

Viscosity::Viscosity(PlasmaDomain &pd): Module(pd) {}

void Viscosity::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "viscosity_output_to_file") m_output_to_file = (rhs[i] == "true");
        else if(lhs[i] == "visc_strength") m_inp_visc_strength = rhs[i];
        else if(lhs[i] == "vars_to_differentiate") m_inp_vars_to_differentiate = rhs[i];
        else if(lhs[i] == "vars_to_update") m_inp_vars_to_update = rhs[i];
        else if(lhs[i] == "visc_opt") m_inp_visc_opt = rhs[i];
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void Viscosity::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if (m_output_to_file) {
    }
}