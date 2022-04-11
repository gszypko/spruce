//sgfilter.hpp
//Header for the Savitzky-Golay Filter Module,
//an implementation of the abstract Module class
//Applies a two-dimensional Savitzky-Golay moving-window
//filter. This results in smoothing of the simulation domain
//according to third-order polynomial fits in each direction.
//Filtering applied repeatedly, after a fixed number of iterations.

#ifndef SGFILTER_HPP
#define SGFILTER_HPP

#include "module.hpp"

class PlasmaDomain;

class SGFilter : public Module {
    public:
        SGFilter(PlasmaDomain &pd);
        void postIterateModule(double dt) override;
        std::string commandLineMessage() const override;
    private:
        void applyFilter(std::vector<Grid> &grids);
        void singleVarSavitzkyGolay(Grid &grid);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        int filter_interval;
};

#endif
