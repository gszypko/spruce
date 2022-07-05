#ifndef EIC_THERMALIZATION_HPP
#define EIC_THERMALIZATION_HPP

#include "module.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "idealmhd2E.hpp"
#include "ideal2F.hpp"
#include <iostream>
#include <functional>
#include <cmath>
#include "constants.hpp"
#include <unordered_map>

class PlasmaDomain;

class EICThermalization : public Module {
    public:
        // *** Construction
        EICThermalization(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        // *** Usage
        void setupModule() override;
        void postIterateModule(double dt) override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
        std::vector<std::string> config_names() const override {return {"eic_output_to_file"};};
    private:
        // *** Members
        std::vector<std::string> m_eqset_grids {"thermal_energy_i","thermal_energy_e","n","temp_e"};
        std::vector<std::string> m_var_names {"a","w_pe","Gam_e","Lam_e","gam_ei","nu_ei","dEdt","dEdEi"};
        enum Vars {a,w_pe,Gam_e,Lam_e,gam_ei,nu_ei,dEdt,dEdEi,num_vars};
        std::vector<Grid> m_vars;
        bool m_output_to_file;
};

#endif