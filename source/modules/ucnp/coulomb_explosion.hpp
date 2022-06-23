#ifndef COULOMB_EXPLOSION_HPP
#define COULOMB_EXPLOSION_HPP

#include "module.hpp"
#include "grid.hpp"
#include "plasmadomain.hpp"
#include "idealmhd.hpp"
#include <iostream>
#include <functional>
#include <cmath>
#include "constants.hpp"

class PlasmaDomain;

class CoulombExplosion : public Module {
    public:
        // *** Construction
        CoulombExplosion(PlasmaDomain &pd);
        // *** Usage
        Grid compute_charge_density(const Grid& r,const Grid& n) const;
        Grid compute_total_charge(const Grid& r,const Grid& rho_c) const;
        void setupModule() override;
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
    private:
        // *** Functions
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        // *** Members
        double m_timescale; // strength of non-neutrality decays exponentially with time
        double m_lengthscale; // length scale of Gaussian non-neutrality function
        double m_strength; // level of non-neutrality at the center
};

#endif