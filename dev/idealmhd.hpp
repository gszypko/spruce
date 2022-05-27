//idealmhd.hpp
//Defines IdealMHD EquationSet

#ifndef IDEALMHD_HPP
#define IDEALMHD_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include <vector>
#include <string>
// Leave BoundaryConditions inside of PlasmaDomain, include way for EquationSet to categorize quantities (as momenta, thermal energies, fields, etc.)
// so that the boundary extrapolations can simply iterate over all of the quantities of a particular type

// d_x,d_y and pos_x,pos_y should remain in the PlasmaDomain, not split out
// Also leave be_x,be_y as part of the plasmadomain, bi_x,bi_y should be part of EquationSet

// Make a getter for the m_grids, by Enum, by integer, or by name string
// Make the Enum and variable names publicly accessible (Enum cannot be defined as "abstract" in the sense that definition is enforced; a redefined enum in a derived class is an entirely separate entity)
// AND make getters for densities, momenta, energies, etc. (for purposes of boundary conditions)
// Make num_variables into a function call

//NOTE: NEED TO HANDLE recomputeTemperature, which refers to a specific variable by definition
//Maybe should instead have a generic "makeVariablesConsistent" function, with recomputeTemperature as a function particular to IdealMHD


class PlasmaDomain;

class IdealMHD: public EquationSet {
    public:
        IdealMHD(PlasmaDomain &pd);
                
        std::vector<int> state_variables() override;
        std::vector<int> densities() override;
        std::vector<int> momenta() override;
        std::vector<int> thermal_energies() override;
        std::vector<int> fields() override;

        void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) override;
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) override;
        Grid computeMagneticEnergyTerm() override;
        void computeConstantTerms() override;
        void recomputeDerivedVariables(std::vector<Grid> &grids) override;
        void recomputeTemperature(std::vector<Grid> &grids); //need to rethink this in the general case
        void catchUnderdensity(std::vector<Grid> &grids) override;
        void recomputeDT() override;
    
    private:
    const std::vector<std::string> m_var_names{"d_x","d_y","pos_x","pos_y","rho","temp","mom_x","mom_y","be_x","be_y","bi_x","bi_y","grav_x","grav_y",
                "press","thermal_energy","kinetic_energy","div_bi","dt","v_x","v_y","n",
                "div_be","b_hat_x","b_hat_y","b_magnitude","v_a"};

};

#endif