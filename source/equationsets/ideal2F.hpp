//idealmhd2F.hpp
//Defines Ideal2F EquationSet

#ifndef IDEAL_2F_HPP
#define IDEAL_2F_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include "constants.hpp"
#include "savitzkygolay.hpp"
#include <vector>
#include <string>

class PlasmaDomain;

class Ideal2F: public EquationSet {
    public:
        Ideal2F(PlasmaDomain &pd);
        Ideal2F() = default;
        std::vector<std::string> config_names() const override 
            {return {"use_sub_cycling","epsilon_courant","smooth_fields","viscosity"};};
        std::vector<std::string> def_var_names() const override{
            return {"i_rho","e_rho","i_mom_x","i_mom_y","e_mom_x","e_mom_y",
                "i_temp","e_temp","bi_x","bi_y","bi_z","E_x","E_y","E_z","grav_x","grav_y",
                "i_n","e_n","i_v_x","i_v_y","e_v_x","e_v_y","j_x","j_y",
                "i_press","e_press","press","i_thermal_energy","e_thermal_energy",
                "rho","rho_c","n","dn","dt",
                "b_x","b_y","b_z","b_mag","b_mag_xy","b_hat_x","b_hat_y",
                "curlE_z","divE","divB","j_x_sg","E_x_sg","E_y_sg","curlE_z_sg","B_z_sg"};
        }
        enum Vars {i_rho,e_rho,i_mom_x,i_mom_y,e_mom_x,e_mom_y,
                i_temp,e_temp,bi_x,bi_y,bi_z,E_x,E_y,E_z,grav_x,grav_y,
                i_n,e_n,i_v_x,i_v_y,e_v_x,e_v_y,j_x,j_y,
                i_press,e_press,press,i_thermal_energy,e_thermal_energy,
                rho,rho_c,n,dn,dt,
                b_x,b_y,b_z,b_mag,b_mag_xy,b_hat_x,b_hat_y,
                curlE_z,divE,divB,j_x_sg,E_x_sg,E_y_sg,curlE_z_sg,B_z_sg};

        std::vector<int> state_variables() override {
            return {i_rho,e_rho,i_mom_x,i_mom_y,e_mom_x,e_mom_y,i_temp,e_temp,bi_x,bi_y,bi_z,E_x,E_y,E_z,grav_x,grav_y};
        }
        std::vector<int> densities() override { return {i_rho,e_rho}; }
        std::vector<std::vector<int>> momenta() override { return {{i_mom_x,i_mom_y},{e_mom_x,e_mom_y}}; }
        std::vector<int> thermal_energies() override { return {i_thermal_energy,e_thermal_energy}; }
        std::vector<int> fields() override { return {}; }

        Grid getDT() override {return m_grids[dt];};

        void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) override;
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) override;
        void populateVariablesFromState(std::vector<Grid> &grids) override;
        void propagateChanges(std::vector<Grid> &grids) override;

    private:
        // private members
        bool m_use_sub_cycling{true};
        double m_epsilon_courant{0.1};
        bool m_verbose{false};

        // private functions
        void recomputeDT();
        void catchNullFieldDirection(std::vector<Grid> &grids);
        void enforceMinimums(std::vector<Grid>& grids);
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids);
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids);
        std::vector<Grid> subcycleMaxwell(const std::vector<Grid>& grids, const std::vector<Grid>& dj, double step);
        void maxwellCurlEqs(const std::vector<Grid>& EM,const std::vector<Grid>& j, std::vector<Grid>& dEM_dt);
        void apply_fixed_curl_bc(Grid& grid) const;
        void smooth_vars(std::vector<Grid> &grids) const;
        void populate_boundary(Grid& grid) const;
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif