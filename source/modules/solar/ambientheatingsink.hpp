//ambientheatingsink.hpp
//Header for the Ambient Heating Sink Module,
//an implementation of the abstract Module class
//Applies a constant, negative volumetric heating rate to a localized region of the domain
//For use in conjunction with the Ambient Heating Module, to reduce ambient heating in localized regions

#ifndef AMBIENTHEATINGSINK_HPP
#define AMBIENTHEATINGSINK_HPP

#include "module.hpp"

class PlasmaDomain;

class AmbientHeatingSink : public Module {
    public:
        AmbientHeatingSink(PlasmaDomain &pd);
        void setupModule() override;
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
    private:
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        double heating_rate{0.0};
        bool exp_mode{false};
        double exp_base_heating_rate;
        double exp_scale_height;
        double center_x; //horizontal position of heating reduction, in cm
        double half_width; //half of the total width of the heating reduction (with parabolic profile falloff), in cm
        double ms_electron_heating_fraction{0.5}; //fraction of direct heating energy given to the electrons, for multispecies analysis (remainder given to ions)
        Grid reduction;     
};

#endif
