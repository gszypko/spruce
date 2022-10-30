#ifndef VISCOSITY_HPP
#define VISCOSITY_HPP

#include "module.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include <vector>
#include <string>

class PlasmaDomain;

class Viscosity : public Module {
    public:
        // *** Construction
        Viscosity(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        // *** Usage
        std::vector<std::string> config_names() const override {return {"viscosity_output_to_file","visc_strength","vars_to_differentiate","vars_to_update"};};
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
    private:
        // *** Members
        std::string m_inp_visc_strength{".01",".01"};
        std::string m_inp_vars_to_differentiate{"v_x","v_y"};
        std::string m_inp_vars_to_update{"mom_x","mom_y"};
        std::string m_inp_visc_opt{"local","local"};
        bool m_output_to_file;
};

#endif