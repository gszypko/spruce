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
    private:
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        double heating_rate{0.0};
        bool exp_mode{false};
        double exp_base_heating_rate;
        double exp_scale_height;
        Grid heating;     
};

#endif
