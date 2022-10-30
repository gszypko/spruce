//idealmhd.hpp
//Defines IdealMHD2E EquationSet

#ifndef IDEALMHD_2E_HPP
#define IDEALMHD_2E_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include <vector>
#include <string>

class PlasmaDomain;

class IdealMHD2E: public EquationSet {
    public:
        IdealMHD2E(PlasmaDomain &pd);
        IdealMHD2E() = default;
        std::vector<std::string> def_var_names() const override{
            return {"rho","i_temp","e_temp","mom_x","mom_y","bi_x","bi_y","grav_x","grav_y",
                "n","i_press","e_press","press","i_thermal_energy","e_thermal_energy","v_x","v_y","kinetic_energy",
                "b_x","b_y","b_mag","b_hat_x","b_hat_y","dt"};
        }
        enum Vars {rho,i_temp,e_temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y,
                n,i_press,e_press,press,i_thermal_energy,e_thermal_energy,v_x,v_y,kinetic_energy,
                b_x,b_y,b_mag,b_hat_x,b_hat_y,dt};

        std::vector<int> state_variables() override {
            return {rho,i_temp,e_temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y};
        }

        std::vector<int> evolved_variables() override {
            return {rho,mom_x,mom_y,i_thermal_energy,e_thermal_energy,bi_x,bi_y};
        }

        std::vector<int> densities() override { return {rho,rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y},{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {i_thermal_energy,e_thermal_energy}; }
        std::vector<int> fields() override { return {bi_x,bi_y}; }

        void computeTimeDerivatives(const std::vector<Grid> &grids) override;
        void applyTimeDerivatives(std::vector<Grid> &grids, double step) override;
        void propagateChanges(std::vector<Grid> &grids) override;
        void viscosity_mask(Grid& grid) const;

    private:
        double m_global_viscosity{0};
        std::string m_viscosity_opt{"velocity"};
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) override;
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) override;
        void recomputeDT() override;
        void catchNullFieldDirection(std::vector<Grid> &grids);
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif