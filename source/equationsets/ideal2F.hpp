//idealmhd2F.hpp
//Defines Ideal2F EquationSet

#ifndef IDEAL_2F_HPP
#define IDEAL_2F_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include <vector>
#include <string>
#include "constants.hpp"

class PlasmaDomain;

class Ideal2F: public EquationSet {
    public:
        Ideal2F(PlasmaDomain &pd);
        Ideal2F() = default;
        std::vector<std::string> def_var_names() const override{
            return {"i_rho","e_rho","i_mom_x","i_mom_y","e_mom_x","e_mom_y",
                "i_temp","e_temp","bi_x","bi_y","bi_z","grav_x","grav_y",
                "i_thermal_energy","e_thermal_energy","E_x","E_y","E_z",
                "i_n","e_n","i_v_x","i_v_y","e_v_x","e_v_y","j_x","j_y",
                "i_press","e_press","press",
                "b_x","b_y","b_z","b_mag","b_mag_xy","b_hat_x","b_hat_y",
                "rho","rho_c","n","dn","divBcond","divEcond","dt"};
        }
        enum Vars {i_rho,e_rho,i_mom_x,i_mom_y,e_mom_x,e_mom_y,
                i_temp,e_temp,bi_x,bi_y,bi_z,grav_x,grav_y,
                i_thermal_energy,e_thermal_energy,E_x,E_y,E_z,
                i_n,e_n,i_v_x,i_v_y,e_v_x,e_v_y,j_x,j_y,
                i_press,e_press,press,
                b_x,b_y,b_z,b_mag,b_mag_xy,b_hat_x,b_hat_y,
                rho,rho_c,n,dn,divBcond,divEcond,dt};

        std::vector<int> state_variables() override {
            return {i_rho,e_rho,i_mom_x,i_mom_y,e_mom_x,e_mom_y,i_temp,e_temp,bi_x,bi_y,bi_z,grav_x,grav_y};
        }
        std::vector<int> densities() override { return {i_rho,e_rho}; }
        std::vector<std::vector<int>> momenta() override { return {{i_mom_x,i_mom_y},{e_mom_x,e_mom_y}}; }
        std::vector<int> thermal_energies() override { return {i_thermal_energy,e_thermal_energy}; }

        Grid getDT() override {return m_grids[dt];};

        void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) override;
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) override;
        void populateVariablesFromState(std::vector<Grid> &grids) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        void recomputeDT();
        void catchNullFieldDirection(std::vector<Grid> &grids);
        void enforceMinimums(std::vector<Grid>& grids);
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids);
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids);
};

#endif