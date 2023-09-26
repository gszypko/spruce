
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
        void iterateModule(double dt) override;
    private:
        double time_scale;
        double safety_factor{1.0};
        Grid diffusivity;
        bool output_to_file = false;
        Grid anomalous_template;
        Grid joule_heating;
        Grid avg_heating;
        double metric_coeff; // = 1.0e49;
        bool metric_smoothing; // = true;
        double smoothing_sigma;
        int kernel_radius;
        int curr_num_subcycles;
        Grid smoothing_kernel;
        void computeDiffusion();
        void computeTemplate(const Grid &b_x, const Grid &b_y);
};

#endif


