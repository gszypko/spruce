//massinjection.hpp
//Header for the Mass Injection Module,
//an implementation of the abstract Module class
//Applies a (Gaussian) localized mass injection rate in the domain
//Location, strength, and time duration specified in module configuration

#ifndef MASSINJECTION_HPP
#define MASSINJECTION_HPP

#include "module.hpp"

class PlasmaDomain;

class MassInjection : public Module {
    public:
        MassInjection(PlasmaDomain &pd);
        void preIterateModule(double dt) override;
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
    private:
        double start_time;
        double duration;
        double max_injection_rate;
        double stddev_x, stddev_y;
        double center_x, center_y;
        Grid injection_template;
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif