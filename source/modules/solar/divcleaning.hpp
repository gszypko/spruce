//divcleaning.hpp
//Header for the Divergence Cleaning Module,
//an implementation of the abstract Module class
//Implements a diffusive divegence cleaning method

#ifndef DIVCLEANING_HPP
#define DIVCLEANING_HPP

#include "module.hpp"

class PlasmaDomain;

class DivCleaning : public Module {
    public:
        DivCleaning(PlasmaDomain &pd);
        void setupModule() override;
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
    private:
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        int computeNumSubcycles(double dt);
        double epsilon{0.1};
        double time_scale{1.0};
        Grid coeff;     
};

#endif
