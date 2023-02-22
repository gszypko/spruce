#ifndef EIC_THERMALIZATION_HPP
#define EIC_THERMALIZATION_HPP

#include "module.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
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
        void computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt);
    private:
        // *** Members
        std::vector<std::string> m_eqset_grids {"i_thermal_energy","e_thermal_energy","rho"};
        std::vector<std::string> m_var_names {"a","w_pe","Gam_e","Lam_e","gam_ei","nu_ei","dEdt","dEdEi"};
        enum Vars {a,w_pe,Gam_e,Lam_e,gam_ei,nu_ei,dEdt,dEdEi,num_vars};
        std::vector<Grid> m_vars;
};

#endif