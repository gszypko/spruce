
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
        void iterateModule(double dt) override;
    private:
        std::string time_integrator{"euler"};
        double safety_factor{1.0};

        std::string template_mode{"flood_fill"};
        double frobenius_metric_coeff{1.0e50};
        double flood_fill_threshold{1.0};
        double flood_fill_current_ramp_length{1.0e-5};
        std::vector<int> null_location;
        double flood_fill_argmin_radius{5.0e9};

        bool metric_smoothing{true}; // = true;
        double smoothing_sigma{3.0};

        std::string resistivity_model{"time_scale"};
        double time_scale{1.0};
        std::vector<double> resistivity_model_params;

        bool gradient_correction{false};

        Grid diffusivity;
        bool output_to_file = false;
        Grid anomalous_template;
        double flood_fill_max_radius{-1.0};
        double flood_fill_min_current{-1.0};
        Grid joule_heating;
        Grid avg_heating;
        int kernel_radius;  
        int curr_num_subcycles;
        Grid smoothing_kernel;
        
        std::vector<int> argminLocalized(std::vector<int> center, const Grid &quantity, const Grid &pos_x, const Grid &pos_y, double radius);
        Grid circularMask(double radius, int i_center, int j_center);
        Grid currentThresholdMask(double threshold, double ramp_length, const Grid &curr_density);
        void computeDiffusion(const Grid &bi_x, const Grid &bi_y, const Grid &bi_z);
        void computeNumSubcycles();
        void computeTemplate(const Grid &b_x, const Grid &b_y, const Grid &curr_density);
        std::vector<Grid> computeTimeDerivatives(const Grid &bi_x, const Grid &bi_y, const Grid &bi_z, const Grid &be_x_lap, const Grid &be_y_lap, const Grid &be_z_lap);
        Grid computeCurrentDensity(const Grid &bi_x, const Grid &bi_y, const Grid &bi_z);
        Grid computeCurrentDensity();
};

#endif


