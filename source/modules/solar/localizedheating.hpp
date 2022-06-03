//localizedheating.hpp
//Header for the Localized Heating Module,
//an implementation of the abstract Module class
//Applies a (Gaussian) localized volumetric heating rate in the domain
//Location, strength, and time duration specified in module configuration

#ifndef LOCALIZEDHEATING_HPP
#define LOCALIZEDHEATING_HPP

#include "module.hpp"

class PlasmaDomain;

class LocalizedHeating : public Module {
    public:
        LocalizedHeating(PlasmaDomain &pd);
        void preIterateModule(double dt) override;
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
    private:
        double start_time;
        double duration;
        double max_heating_rate;
        double stddev_x, stddev_y;
        double center_x, center_y;
        Grid heating_template;
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif