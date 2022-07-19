//idealmhd.hpp
//Defines IdealMHD EquationSet

#ifndef IDEALMHD_HPP
#define IDEALMHD_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include <vector>
#include <string>

class PlasmaDomain;

class IdealMHD: public EquationSet {
    public:
        IdealMHD(PlasmaDomain &pd);
        IdealMHD() = default;
        std::vector<std::string> def_var_names() const override{
            return {"rho","temp","mom_x","mom_y","bi_x","bi_y","grav_x","grav_y",
                "n","press","thermal_energy","v_x","v_y","kinetic_energy",
                "b_x","b_y","b_mag","b_hat_x","b_hat_y","dt"};
        }
        enum Vars {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y,
                n,press,thermal_energy,v_x,v_y,kinetic_energy,
                b_x,b_y,b_mag,b_hat_x,b_hat_y,dt};

        std::vector<int> state_variables() override {
            return {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y};
        }
        std::vector<int> densities() override { return {rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {thermal_energy}; }

        void populateVariablesFromState(std::vector<Grid> &grids) override;
        Grid getDT() override;
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) override;
        void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        void recomputeDT();
        void enforceMinimums(std::vector<Grid> &grids);
        void computeConstantTerms(std::vector<Grid> &grids);
        void recomputeDerivedVariables(std::vector<Grid> &grids);
        void recomputeTemperature(std::vector<Grid> &grids);
        void recomputeThermalEnergy(std::vector<Grid> &grids);
        void recomputeKineticEnergy(std::vector<Grid> &grids);
        void recomputeNumberDensity(std::vector<Grid> &grids);
        void recomputeMagneticFields(std::vector<Grid> &grids);
        std::vector<Grid> computeTimeDerivativesCharacteristicBoundary(const std::vector<Grid> &grids, bool x_bound_1, bool x_bound_2, bool y_bound_1, bool y_bound_2);
        std::vector<std::vector<Grid>> singleBoundaryTermsMOC(const std::vector<Grid> &grids, int boundary_index, bool boundary_lower);
};

#endif