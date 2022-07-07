//idealmhd2F.hpp
//Defines Ideal2F EquationSet

#ifndef IDEAL_2F_HPP
#define IDEAL_2F_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include <vector>
#include <string>

class PlasmaDomain;

class Ideal2F: public EquationSet {
    public:
        Ideal2F(PlasmaDomain &pd);
        Ideal2F() = default;
        std::vector<std::string> def_var_names() const override{
            return {"i_rho","e_rho","i_temp","e_temp","i_mom_x","i_mom_y","e_mom_x","e_mom_y",
                "bi_x","bi_y","grav_x","grav_y","i_n","e_n","i_press","e_press","press","i_thermal_energy","e_thermal_energy",
                "i_v_x","i_v_y","e_v_x","e_v_y","i_kinetic_energy","e_kinetic_energy",
                "E_x","E_y","E_z","b_x","b_y","b_mag","b_hat_x","b_hat_y","dt"};
        }
        enum Vars {i_rho,e_rho,i_temp,e_temp,i_mom_x,i_mom_y,e_mom_x,e_mom_y,
                bi_x,bi_y,grav_x,grav_y,i_press,e_press,press,i_thermal_energy,e_thermal_energy,
                i_n,e_n,i_v_x,i_v_y,e_v_x,e_v_y,i_kinetic_energy,e_kinetic_energy,
                E_x,E_y,E_z,b_x,b_y,b_mag,b_hat_x,b_hat_y,dt};

        std::vector<int> state_variables() override {
            return {i_rho,e_rho,i_temp,e_temp,i_mom_x,i_mom_y,e_mom_x,e_mom_y,bi_x,bi_y,grav_x,grav_y};
        }
        std::vector<int> densities() override { return {i_rho,e_rho}; }
        std::vector<std::vector<int>> momenta() override { return {{i_mom_x,i_mom_y},{e_mom_x,e_mom_y}}; }
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
        void recomputeElectricFields(std::vector<Grid> &grids);
        void catchUnderdensity(std::vector<Grid> &grids);
};

#endif