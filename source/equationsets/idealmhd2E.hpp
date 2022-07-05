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
            return {"rho","temp_i","temp_e","mom_x","mom_y","bi_x","bi_y","grav_x","grav_y",
                "press_i","press_e","press","thermal_energy_i","thermal_energy_e","kinetic_energy",
                "div_bi","dt","v_x","v_y","n","div_be","b_hat_x","b_hat_y","b_magnitude"};
        }
        enum Vars {rho,temp_i,temp_e,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y,
                press_i,press_e,press,thermal_energy_i,thermal_energy_e,kinetic_energy,
                div_bi,dt,v_x,v_y,n,div_be,b_hat_x,b_hat_y,b_magnitude};

        std::vector<int> state_variables() override {
            return {rho,temp_i,temp_e,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y};
        }
        std::vector<int> densities() override { return {rho,rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y},{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {thermal_energy_i,thermal_energy_e}; }

        void populateVariablesFromState(std::vector<Grid> &grids) override;
        Grid getDT() override;
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) override;
        void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        void recomputeDT();
        void computeConstantTerms(std::vector<Grid> &grids);
        void recomputeDerivedVariables(std::vector<Grid> &grids);
        void recomputeTemperature(std::vector<Grid> &grids); //need to rethink this in the general case
        void catchUnderdensity(std::vector<Grid> &grids);
};

#endif