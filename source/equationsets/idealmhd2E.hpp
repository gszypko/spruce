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
        std::vector<int> densities() override { return {rho,rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y},{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {i_thermal_energy,e_thermal_energy}; }

        void populateVariablesFromState(std::vector<Grid> &grids) override;
        Grid getDT() override;
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) override;
        void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        void recomputeDT();
        void computeConstantTerms(std::vector<Grid> &grids);
        void recomputeDerivedVariables(std::vector<Grid> &grids);
        void recomputeTemperature(std::vector<Grid> &grids);
        void recomputeThermalEnergy(std::vector<Grid> &grids);
        void recomputeKineticEnergy(std::vector<Grid> &grids);
        void recomputeNumberDensity(std::vector<Grid> &grids);
        void recomputeMagneticFields(std::vector<Grid> &grids);
        void enforceMinimums(std::vector<Grid>& grids);
};

#endif