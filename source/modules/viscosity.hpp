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
        bool is_momentum(std::string name) const; // check whether variable name is a momentum variable
        bool is_velocity(std::string name) const; // check whether variable name is a velocity variable
        bool is_thermal_energy(std::string name) const; // check whether variable name is a velocity variable
        bool is_temperature(std::string name) const; // check whether variable name is a velocity variable
        // *** Construction
        Viscosity(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override; // handles configs during construction
        void setupModule(); // to be called in PlasmaDomain constructor - does initial setup of module following config parsing
        // *** Usage
        Grid constructSingleViscosityGrid(const std::vector<Grid>& grids, int term_index);
        std::vector<Grid> constructViscosityGrids(const std::vector<Grid>& grids);
        Grid getBoundaryViscosity(double strength,double length) const;
        void computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt);
        void iterateModule(double dt) override;
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
        std::vector<Grid> m_grids_strength;
        std::vector<Grid> m_grids_lap;
        std::vector<Grid> m_grids_dqdt;
        std::vector<std::string> m_dt_names; // <name>_dt
        std::vector<std::string> m_strength_names; // <name>_str
        std::vector<std::string> m_lap_names; // <name>_lap
        std::vector<std::string> m_dqdt_names; // <name>_dqdt

        bool m_output_visc;
        bool m_output_lap;
        bool m_output_strength;
        bool m_output_timescale;

        std::string m_hv_time_integrator; // time integrator used for hyperviscosity (strength > 1.0), one of {euler, rk2, rk4}
        double m_hv_epsilon{1.0}; // epsilon used for subcycling in hyperviscous (strength > 1.0) evolution

        std::vector<std::string> m_momenta; // possible momentum density variables
        std::vector<std::string> m_velocities; // possible velocity variables
        std::vector<std::string> m_thermal_energies; // possible thermal energy variables
        std::vector<std::string> m_temperatures; // possible temperature variables
};

#endif