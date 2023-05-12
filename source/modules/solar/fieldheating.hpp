//fieldheating.hpp
//Header for the Field Heating Module,
//an implementation of the abstract Module class
//Applies volumetric heating to the plasma proportional
//to the magnitude of the local magnetic field

#ifndef FIELDHEATING_HPP
#define FIELDHEATING_HPP

#include "module.hpp"

class FieldHeating : public Module {
    public:
        FieldHeating(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
        void setupModule() override;
    private:
        // double rate;
        double coeff = 0.0, current_pow = 0.0, b_pow = 0.0, n_pow = 0.0, roc_pow = 0.0;
        // std::string mode;
        bool output_to_file;
        // bool current_mode; //when true, heating is calculated relative to current density J instead of magnetic pressure
        Grid heating;
        void computeHeating();
};

#endif


