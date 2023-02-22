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
                "b_x","b_y","b_mag","b_hat_x","b_hat_y","dt","dPdx","a_x","e_temp_g"};
        }
        enum Vars {rho,i_temp,e_temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y,
                n,i_press,e_press,press,i_thermal_energy,e_thermal_energy,v_x,v_y,kinetic_energy,
                b_x,b_y,b_mag,b_hat_x,b_hat_y,dt};

        std::vector<int> state_variables() const override {
            return {rho,i_temp,e_temp,mom_x,mom_y,bi_x,bi_y,grav_x,grav_y};
        }

        std::vector<int> evolved_variables() const override {
            return {rho,mom_x,mom_y,i_thermal_energy,e_thermal_energy,bi_x,bi_y};
        }

        std::vector<std::string> species() const override {return {"i","e"};}
        std::vector<int> densities() const override { return {rho,rho}; }
        std::vector<int> number_densities() const override { return {n,n}; }
        std::vector<std::vector<int>> momenta() const override { return {{mom_x,mom_y},{mom_x,mom_y}}; }
        std::vector<std::vector<int>> velocities() const override { return {{v_x,v_y},{v_x,v_y}}; }
        std::vector<int> thermal_energies() const override { return {i_thermal_energy,e_thermal_energy}; }
        std::vector<int> pressures() const override { return {i_press,e_press}; }
        std::vector<int> temperatures() const override { return {i_temp,e_temp}; }
        std::vector<int> fields() const override { return {bi_x,bi_y}; }
        std::vector<int> timescale() const override {return {dt,dt}; }

    private:
        double m_global_viscosity{0};
        std::string m_viscosity_opt{"velocity"};

        std::vector<Grid> computeTimeDerivativesDerived(const std::vector<Grid> &grids) const override;
        void enforceMinimums(std::vector<Grid>& grids) const override;
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) const override;
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) const override;
        void recomputeDT(std::vector<Grid>& grids) const override;
        void catchNullFieldDirection(std::vector<Grid> &grids) const;
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif