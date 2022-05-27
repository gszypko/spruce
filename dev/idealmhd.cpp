#include "idealmhd.hpp"

IdealMHD::IdealMHD(PlasmaDomain &pd): EquationSet(pd,def_var_names()) {}

std::vector<int> IdealMHD::state_variables(){ return {0}; }
std::vector<int> IdealMHD::densities(){ return{0};  }
std::vector<int> IdealMHD::momenta(){ return{0};  }
std::vector<int> IdealMHD::thermal_energies(){ return{0};  }
std::vector<int> IdealMHD::fields(){ return{0};  }

void IdealMHD::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step){ return; }
std::vector<Grid> IdealMHD::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff){ return {Grid::Zero(1,1)}; }
Grid IdealMHD::computeMagneticEnergyTerm(){ return Grid::Zero(1,1); }
void IdealMHD::computeConstantTerms(){ return; }
void IdealMHD::recomputeDerivedVariables(std::vector<Grid> &grids){ return; }
void IdealMHD::recomputeTemperature(std::vector<Grid> &grids){ return; } //need to rethink this in the general case
void IdealMHD::catchUnderdensity(std::vector<Grid> &grids){ return; }
void IdealMHD::recomputeDT(){ return; }