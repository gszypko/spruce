
#ifndef ANOMALOUSRESISTIVITY_HPP
#define ANOMALOUSRESISTIVITY_HPP

#include "module.hpp"

class AnomalousResistivity : public Module {
    public:
        AnomalousResistivity(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void setupModule() override;
        std::string commandLineMessage() const override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
        void computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt) override;
    private:
        double time_scale;
        Grid diffusivity;
        bool output_to_file = false;
        Grid anomalous_template;
        double metric_coeff; // = 1.0e49;
        bool metric_smoothing; // = true;
        double smoothing_sigma;
        int kernel_radius;
        Grid smoothing_kernel;
        void computeDiffusion(const std::vector<Grid> &grids);
        void computeTemplate(const std::vector<Grid> &grids);
};

#endif


