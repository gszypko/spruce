
#ifndef MOMENTUMINJECTION_HPP
#define MOMENTUMINJECTION_HPP

#include "module.hpp"
#include <vector>

class PlasmaDomain;

class MomentumInjection : public Module {
    public:
        MomentumInjection(PlasmaDomain &pd);
        void postIterateModule(double dt) override;
        void setupModule() override;
        std::string commandLineMessage() const override;
    private:
        bool oscillatory{false};
        double oscillation_period;
        double start_time;
        double duration;
        double max_accel; //maximum acceleration to apply, reached at center of gaussian injection template
        double stddev_x, stddev_y; //std devs of injection template
        double center_x, center_y; //center coords of injection template, in grid coords
        double template_angle{0.0}; //angle, in degrees, to rotate the template clockwise (does NOT affect acceleration direction)
        std::vector<double> dir; //components of a vector describing the direction of acceleration; will be normalized automatically
        std::vector<Grid> accel_template;
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
};

#endif