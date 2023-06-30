//radiativelosses.hpp
//Header for the Radiative Losses Module,
//an implementation of the abstract Module class
//Applies a piecewise power-law approximation
//of optically-thin radiation in the solar corona

#ifndef RADIATIVELOSSES_HPP
#define RADIATIVELOSSES_HPP

#include "module.hpp"

class RadiativeLosses : public Module {
    public:
        RadiativeLosses(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void preIterateModule(double dt) override;
        void iterateModule(double dt) override;
        std::string commandLineMessage() const override;
        void setupModule() override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
    private:
        double cutoff_ramp;
        double cutoff_temp;
        double epsilon;
        bool output_to_file = false;
        Grid avg_losses;
        int curr_num_subcycles;
        std::string time_integrator;
        bool inactive_mode = false; //when true, module will compute the change in thermal energy (allowing for output), but WON'T apply it to the simulation
        bool prevent_subcycling = false; //when true, loss rates won't be allowed to get high enough to dominate the fluid time scales
        Grid computeLosses(const Grid &temp, const Grid &n) const;
        Grid computeLosses() const;
        int numberSubcycles(double dt);
};

#endif


