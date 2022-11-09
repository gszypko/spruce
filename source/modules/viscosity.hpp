#ifndef VISCOSITY_HPP
#define VISCOSITY_HPP

#include "module.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include <vector>
#include <string>
#include "utils.hpp"
#include "equationset.hpp"
#include "constants.hpp"

class PlasmaDomain;

class Viscosity : public Module {
    public:
        // *** Variable Types
        std::vector<std::string> m_momenta {"mom_x","mom_y","i_mom_x","i_mom_y","e_mom_x","e_mom_y"}; // possible momentum density variables across all equation sets
        std::vector<std::string> m_velocities {"v_x","v_y","i_v_x","i_v_y","e_v_x","e_v_y"}; // possible velocity variables across all equation sets
        bool is_momentum(std::string name) const; // check whether variable name is a momentum variable
        bool is_velocity(std::string name) const; // check whether variable name is a velocity variable
        // *** Construction
        Viscosity(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override; // handles configs during construction
        void setupModule(); // to be called in PlasmaDomain constructor - does initial setup of module following config parsing
        // *** Usage
        void constructViscosityGrids(const std::vector<Grid>& grids);
        Grid getBoundaryViscosity(double strength,double length) const;
        void computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt);
        std::vector<std::string> config_names() const override {return {"visc_output_to_file","visc_strength","visc_vars_diff","visc_vars_update"};};
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
    private:
        // *** Members
        std::string m_inp_visc_opt; // global, local, or boundary
        std::string m_inp_strength; // strength coefficient for viscosity term
        std::string m_inp_vars_to_diff; // variables to differentiate to compute viscosity term
        std::string m_inp_vars_to_evol; // variables that are evolved with viscosity term
        std::string m_inp_length; // length scale for decay of boundary viscosity - does not apply for other viscosity types, but must be specified anyway
        std::string m_inp_species; // either e = electron, i = ion, or f = fields
        int m_num_terms;

        std::vector<std::string> m_visc_opt;
        std::vector<double> m_strength;
        std::vector<std::string> m_vars_to_diff;
        std::vector<std::string> m_vars_to_evol;
        std::vector<double> m_length;
        std::vector<std::string> m_species;

        std::vector<Grid> m_grids_dt;
        std::vector<Grid> m_grids_lap;
        std::vector<Grid> m_grids_strength;
        std::vector<std::string> m_dt_names;
        std::vector<std::string> m_lap_names;
        std::vector<std::string> m_strength_names;

        bool m_output_visc;
        bool m_output_lap;
        bool m_output_strength;
};

#endif