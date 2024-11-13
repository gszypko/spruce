//viscousheating.hpp
//Header for the Viscous Heating Module,
//an implementation of the abstract Module class
//Applies viscous heating rate as described by
//Braginskii (1965), and Wang et al. 2024

#ifndef VISCOUSHEATING_HPP
#define VISCOUSHEATING_HPP

#include "module.hpp"

class ViscousHeating : public Module {
    public:
        ViscousHeating(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void preIterateModule(double dt) override;
        void iterateModule(double dt) override;
        std::string commandLineMessage() const override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
        void setupModule() override;
    private:
        double coeff = 0.0;
        bool output_to_file = false;
        bool inactive_mode = false; //when true, module will compute the change in thermal energy (allowing for output), but WON'T apply it to the simulation
        Grid heating;
        void computeHeating();
};

#endif


