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
                "b_x","b_y","b_mag","b_hat_x","b_hat_y","dt","dPdx","dPdy"};
        }
        enum Vars {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y,
                n,press,thermal_energy,v_x,v_y,kinetic_energy,
                b_x,b_y,b_mag,b_hat_x,b_hat_y,dt,dPdx,dPdy};

        std::vector<int> state_variables() const override {
            return {rho,temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y};
        }

        std::vector<std::string> evolved_var_names() const override {
            return {"rho","mom_x","mom_y","thermal_energy","bi_x","bi_y"};
        }

        std::vector<int> densities() override { return {rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {thermal_energy}; }
        std::vector<int> fields() override { return {bi_x,bi_y}; }

        void computeTimeDerivativesDerived(const std::vector<Grid> &grids, std::vector<Grid> &grids_dt) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        double m_global_viscosity{0};
        std::string m_viscosity_opt{"velocity"};
        double m_guide_field{};
        void setupEquationSetDerived() override;
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) override;
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) override;
        void catchNullFieldDirection(std::vector<Grid> &grids);
        void recomputeDT() override;
        std::vector<Grid> computeTimeDerivativesCharacteristicBoundary(const std::vector<Grid> &grids, bool x_bound_1, bool x_bound_2, bool y_bound_1, bool y_bound_2, double visc_coeff);
        std::vector<std::vector<Grid>> singleBoundaryTermsMOC(const std::vector<Grid> &grids, int boundary_index, bool boundary_lower, double visc_coeff);
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif