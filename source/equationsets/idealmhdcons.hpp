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

        std::vector<int> evolved_variables() override {
            return {rho,mom_x,mom_y,energy,bi_x,bi_y};
        }

        std::vector<int> densities() override { return {rho}; }
        std::vector<std::vector<int>> momenta() override { return {{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() override { return {thermal_energy}; }
        std::vector<int> fields() override { return {bi_x,bi_y}; }

        void computeTimeDerivatives(const std::vector<Grid> &grids) override;
        void applyTimeDerivatives(std::vector<Grid> &grids, double step) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        double m_global_viscosity{0};
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) override;
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) override;
        void recomputeDT() override;
        void catchNullFieldDirection(std::vector<Grid> &grids);
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif