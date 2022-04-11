//radiativelosses.hpp
//Header for the Radiative Losses Module,
//an implementation of the abstract Module class
//Applies a piecewise power-law approximation
//of optically-thin radiation in the solar corona

#ifndef RADIATIVELOSSES_HPP
#define RADIATIVELOSSES_HPP

#include "module.hpp"
#include "plasmadomain.hpp"

class RadiativeLosses : public Module {
    public:
        RadiativeLosses(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void preIterateModule(double dt) override;
        void iterateModule(double dt) override;
        std::string commandLineMessage() const override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const override;
    private:
        double cutoff_ramp;
        double cutoff_temp;
        double epsilon;
        bool output_to_file;
        Grid next_losses, avg_losses;
        int curr_num_subcycles;
        void computeLosses();
        int numberSubcycles(double dt);
};

#endif


