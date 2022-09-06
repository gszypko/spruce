//idealmhdcons.hpp
//Defines idealmhdcons EquationSet

#ifndef IDEALMHDCONS_HPP
#define IDEALMHDCONS_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include <vector>
#include <string>

class PlasmaDomain;

class IdealMHDCons: public EquationSet {
    public:
        IdealMHDCons(PlasmaDomain &pd);
        IdealMHDCons() = default;
        std::vector<std::string> def_var_names() const override{
            return {"rho","temp","mom_x","mom_y","bi_x","bi_y","grav_x","grav_y",
                "n","v_x","v_y","kinetic_energy","b_x","b_y","b_mag","mag_energy","b_hat_x","b_hat_y",
                "thermal_energy","press","press_tot","energy","dt"};
        }
        enum Vars {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y,
                n,v_x,v_y,kinetic_energy,b_x,b_y,b_mag,mag_energy,b_hat_x,b_hat_y,
                thermal_energy,press,press_tot,energy,dt};

        std::vector<int> state_variables() override {
            return {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y};
        }
        std::vector<int> densities() override { return {rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {thermal_energy}; }
        std::vector<int> fields() override { return {bi_x,bi_y}; }

        Grid getDT() override;
        void populateVariablesFromState(std::vector<Grid> &grids) override;
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids) override;
        void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        double m_global_viscosity{0};
        void enforceMinimums(std::vector<Grid>& grids);
        void recomputeDT();
        void computeConstantTerms(std::vector<Grid> &grids);
        void recomputeDerivedVariables(std::vector<Grid> &grids);
        void recomputeNumberDensity(std::vector<Grid> &grids);
        void recomputeKineticEnergy(std::vector<Grid> &grids);
        void recomputeMagneticEnergy(std::vector<Grid> &grids);
        void recomputeThermalEnergy(std::vector<Grid> &grids);
        void recomputeTemperature(std::vector<Grid> &grids);
        void recomputeThermalEnergyFromTemp(std::vector<Grid> &grids);
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif