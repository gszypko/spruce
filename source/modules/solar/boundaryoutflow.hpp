
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
        Grid accel_template;
        std::string boundary{"y_bound_2"};
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        Grid constructBoundaryAccel(double strength,double length) const;
};

#endif