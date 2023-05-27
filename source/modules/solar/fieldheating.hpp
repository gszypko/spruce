//fieldheating.hpp
//Header for the Field Heating Module,
//an implementation of the abstract Module class
//Applies volumetric heating to the plasma according to
//a user-defined power law heuristic based on magnetic
//parameters and number density

#ifndef FIELDHEATING_HPP
#define FIELDHEATING_HPP

#include "module.hpp"

class FieldHeating : public Module {
    public:
        FieldHeating(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void preIterateModule(double dt) override;
        void iterateModule(double dt) override;
        std::string commandLineMessage() const override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
        void setupModule() override;
    private:
        double coeff = 0.0, current_pow = 0.0, b_pow = 0.0, n_pow = 0.0, roc_pow = 0.0;
        bool output_to_file;
        bool inactive_mode = false; //when true, module will compute the change in thermal energy (allowing for output), but WON'T apply it to the simulation
        Grid heating;
        void computeHeating();
};

#endif


