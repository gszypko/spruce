
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
        double diffusivity;
        Grid anomalous_template;
        void computeTemplate();
        void computeDiffusion();
};

#endif


