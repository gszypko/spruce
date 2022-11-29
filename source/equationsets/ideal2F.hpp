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
        void setupEquationSetDerived() override; 
        std::vector<std::string> config_names() const override 
            {return {"use_sub_cycling","epsilon_courant","smooth_fields","viscosity"};};
        std::vector<std::string> def_var_names() const override{
            return {"i_rho","e_rho","i_mom_x","i_mom_y","e_mom_x","e_mom_y",
                "i_temp","e_temp","bi_x","bi_y","bi_z","E_x","E_y","E_z","grav_x","grav_y",
                "i_n","e_n","i_v_x","i_v_y","e_v_x","e_v_y","j_x","j_y",
                "i_press","e_press","press","i_thermal_energy","e_thermal_energy",
                "rho","rho_c","n","dn","dt","dt_i",
                "b_x","b_y","b_z","b_mag","b_mag_xy","b_hat_x","b_hat_y",
                "curlE_z","divE","divB","lapEx","E_ideal","divEcond"};
        }
        enum Vars {i_rho,e_rho,i_mom_x,i_mom_y,e_mom_x,e_mom_y,
                i_temp,e_temp,bi_x,bi_y,bi_z,E_x,E_y,E_z,grav_x,grav_y,
                i_n,e_n,i_v_x,i_v_y,e_v_x,e_v_y,j_x,j_y,
                i_press,e_press,press,i_thermal_energy,e_thermal_energy,
                rho,rho_c,n,dn,dt,dt_i,
                b_x,b_y,b_z,b_mag,b_mag_xy,b_hat_x,b_hat_y,
                curlE_z,divE,divB,lapEx,E_ideal,divEcond};

        std::vector<int> state_variables() const override {
            return {i_rho,e_rho,i_mom_x,i_mom_y,e_mom_x,e_mom_y,i_temp,e_temp,bi_x,bi_y,bi_z,E_x,E_y,E_z,grav_x,grav_y};
        }

        std::vector<int> evolved_variables() const override {
            return {i_rho,e_rho,i_mom_x,i_mom_y,e_mom_x,e_mom_y,i_thermal_energy,e_thermal_energy,E_x,E_y,E_z,bi_x,bi_y,bi_z};
        }

        std::vector<int> densities() const override { return {i_rho,e_rho}; }
        std::vector<std::vector<int>> momenta() const override { return {{i_mom_x,i_mom_y},{e_mom_x,e_mom_y}}; }
        std::vector<int> thermal_energies() const override { return {i_thermal_energy,e_thermal_energy}; }
        std::vector<int> fields() const override { return {}; }
        std::vector<int> timescale() const override {return {dt_i,dt}; }

    private:
        // private members
        bool m_use_sub_cycling{true};
        double m_epsilon_courant{0.1};
        bool m_verbose{false};
        std::string m_viscosity_opt{"momentum"};
        bool m_remove_curl_terms{false};
        double m_global_viscosity{0};

        // private functions
        std::vector<Grid> computeTimeDerivativesDerived(const std::vector<Grid> &grids) const override;
        void recomputeDT(std::vector<Grid>& grids) const override;
        Grid ionTimescale() const;
        void catchNullFieldDirection(std::vector<Grid> &grids) const;
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) const override;
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) const override;
        std::vector<Grid> subcycleMaxwell(const std::vector<Grid>& grids, const std::vector<Grid>& dj, double step) const;
        void maxwellCurlEqs(const std::vector<Grid>& EM,const std::vector<Grid>& j, std::vector<Grid>& EM_laplacian, std::vector<Grid>& dEM_dt) const;
        void apply_fixed_curl_bc(Grid& grid) const;
        void smooth_vars(std::vector<Grid> &grids) const;
        void populate_boundary(Grid& grid) const;
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif