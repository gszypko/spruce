#ifndef EIC_THERMALIZATION_HPP
#define EIC_THERMALIZATION_HPP

#include "module.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "idealmhd2E.hpp"
#include <iostream>
#include <functional>
#include <cmath>
#include "constants.hpp"
#include <unordered_map>

class PlasmaDomain;

class EICThermalization : public Module {
    public:
        // *** Construction
        EICThermalization(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        // *** Usage
        void setupModule() override;
        void postIterateModule(double dt) override;
        Grid compute_charge_density(const Grid& r,const Grid& n) const;
        Grid compute_total_charge(const Grid& r,const Grid& rho_c) const;
        std::string commandLineMessage() const override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
    private:
        // *** Members
        double m_timescale; // strength of non-neutrality decays exponentially with time
        double m_lengthscale; // length scale of Gaussian non-neutrality function
        double m_strength; // level of non-neutrality at the center
        bool output_to_file{};
        std::vector<std::string> m_required_grids {"press","n","mom_x","mom_y"};
        std::unordered_map<std::string,int> m_grid_ind;
        Grid m_Fcx, m_Fcy, m_dPx, m_dPy;
};

#endif