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
                "press","thermal_energy","kinetic_energy","div_bi","dt","v_x","v_y","n",
                "div_be","b_hat_x","b_hat_y","b_magnitude","dPx","dPy","Fcx","Fcy"};
        }
        enum Vars {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y,
                press,thermal_energy,kinetic_energy,div_bi,dt,v_x,v_y,n,
                div_be,b_hat_x,b_hat_y,b_magnitude,dPx,dPy,Fcx,Fcy};

        std::vector<int> state_variables() override {
            return {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y};
        }
        std::vector<int> densities() override { return {rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {thermal_energy}; }
        std::vector<std::vector<int>> fields() override { return {{bi_x,bi_y}}; }

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
        std::vector<Grid> computeTimeDerivativesCharacteristicBoundary(const std::vector<Grid> &grids, bool x_bound_1, bool x_bound_2, bool y_bound_1, bool y_bound_2);
        std::vector<Grid> singleBoundaryTermsMOC(const std::vector<Grid> &grids, int boundary_index, bool boundary_lower);

};

#endif