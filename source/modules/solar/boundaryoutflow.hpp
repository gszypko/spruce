
#ifndef BOUNDARYOUTFLOW_HPP
#define BOUNDARYOUTFLOW_HPP

#include "module.hpp"
#include <vector>

class PlasmaDomain;

class BoundaryOutflow : public Module {
    public:
        BoundaryOutflow(PlasmaDomain &pd);
        void postIterateModule(double dt) override;
        void setupModule() override;
        std::string commandLineMessage() const override;
    private:
        double max_accel; //maximum acceleration to apply, reached at center of gaussian injection template
        double falloff_length; //scale height of exponential falloff from boundary
        double feather_length{0.0}; //scale length to smooth the edges of the outflow acceleration to zero; ramps over about two times this length
        double mean_outflow{0.0};
        double dynamic_time{1.0};
        double dynamic_target_speed{0.0};
        double curr_accel{0.0};
        Grid accel_template;
        std::string boundary{"y_bound_2"};
        std::string falloff_shape{"exp"};
        bool field_aligned_mode{false};
        bool dynamic_mode{false};
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        Grid constructBoundaryAccel(double strength,double length) const;
        double computeMeanOutflow() const;
};

#endif