//ambientheating.hpp
//Header for the Ambient Heating Module,
//an implementation of the abstract Module class
//Applies a constant, ambient volumetric heating rate to the entire domain

#ifndef AMBIENTHEATING_HPP
#define AMBIENTHEATING_HPP

#include "module.hpp"

class PlasmaDomain;

class AmbientHeating : public Module {
    public:
        AmbientHeating(PlasmaDomain &pd);
        void setupModule() override;
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
        Grid radiativeLosses();
    private:
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        double heating_rate;
        bool rad_mirror_mode = false; //when true, heating is taken as the radiative loss rate for the initial state, constant in time
        double rad_mirror_factor; //coefficient to multiply radiative loss by, when applying as heating
        double cutoff_ramp, cutoff_temp; //radiative losses cutoff for low temps; same as in radiativelosses.hpp
        Grid heating;     
};

#endif
