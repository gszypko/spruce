//physicalviscosity.hpp
//Header for the Physical Viscosity Module,
//an implementation of the abstract Module class
//Applies viscosity and viscous heating as described by
//Braginskii (1965), and Wang et al. 2024 (only eta_0)

#ifndef PHYSICALVISCOSITY_HPP
#define PHYSICALVISCOSITY_HPP

#include "module.hpp"

class PhysicalViscosity : public Module {
    public:
        PhysicalViscosity(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void preIterateModule(double dt) override;
        void iterateModule(double dt) override;
        std::string commandLineMessage() const override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
        void setupModule() override;
    private:
        double coeff = 0.0;
        Grid coeff_grid;
        double ramp_length = 0.0; //physical length of ramp-down from max coefficient to zero, from interior toward boundary
        double buffer_length = 0.0; //physical thickness of exterior border that's set at zero coefficient
        double epsilon = 1.0;
        bool heating_on = true;
        bool force_on = true;
        int num_subcycles;
        bool output_to_file = false;
        bool inactive_mode = false; //when true, module will compute the change in thermal energy (allowing for output), but WON'T apply it to the simulation
        bool gradient_correction = false;
        std::string time_integrator;
        Grid avg_heating;
        std::vector<Grid> avg_force; //average over single timstep, across subcycles
        // Grid computeHeating();
        Grid computeHeating(const std::vector<Grid>& v, const std::vector<Grid>& b_hat, const Grid& temperature);
        int computeViscousSubcycles(double dt);
        std::vector<Grid> computeViscousForce(const std::vector<Grid>& v, const std::vector<Grid>& b_hat, const Grid& temperature);
        Grid constructCoefficientGrid(double strength,double ramp_length,double buffer_length) const;
};

#endif


