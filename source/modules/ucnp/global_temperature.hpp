#ifndef GLOBAL_TEMPERATURE_HPP
#define GLOBAL_TEMPERATURE_HPP

#include "module.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include <iostream>
#include <functional>
#include <cmath>
#include "constants.hpp"
#include <unordered_map>
#include "utils.hpp"

class PlasmaDomain;

class GlobalTemperature : public Module {
    public:
        GlobalTemperature(PlasmaDomain &pd);
        std::vector<std::string> config_names() const override {return {"gt_species","gt_strength","gt_use_diffusion","gt_use_global_temp"};};
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        void setupModule() override;
        void preRecomputeDerivedModule(std::vector<Grid>& grids) const override;
        void postIterateModule(double dt) override;
    private:
        bool m_use_diffusion{false};
        bool m_use_global_temp{false};
        std::vector<std::string> m_species;
        std::vector<int> m_species_ind; // m_species_ind[i] holds the equation set species index for m_species_name[i]
        std::vector<double> m_species_mass; // m_species_mass[i] holds the mass for m_species_name[i]
        double m_strength; // strength factor for diffusion term
        Grid m_dr; // spatial grid for diffusion term
        Grid m_coeff; // diffusion coefficient
        double m_dt_diff; // diffusion timescale
        double m_dt_rk; // rk timestep for diffusion term
        int m_num_rk_steps;
};

#endif