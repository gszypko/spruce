//equationset.hpp
//Defines abstract base class EquationSet
//All specific EquationSets must be defined/implemented as
//a derived class of EquationSet

#ifndef EQUATIONSET_HPP
#define EQUATIONSET_HPP

#include "grid.hpp"
#include <vector>
#include <string>
#include <cassert>
#include <unordered_map>
#include <memory>
// Functions that need to be split out into EquationSet
// evolution.cpp
// void PlasmaDomain::applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step)
// std::vector<Grid> PlasmaDomain::computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff)
// Grid PlasmaDomain::computeMagneticEnergyTerm()
// void PlasmaDomain::computeConstantTerms()
// void PlasmaDomain::recomputeDerivedVariables(std::vector<Grid> &grids)
// void PlasmaDomain::recomputeTemperature(std::vector<Grid> &grids)
// void PlasmaDomain::catchUnderdensity(std::vector<Grid> &grids)
// void PlasmaDomain::recomputeDT()


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

class EquationSet {
    public:
        EquationSet(PlasmaDomain &pd, std::vector<std::string> var_names);
        virtual ~EquationSet() {}
        static std::unique_ptr<EquationSet> spawnEquationSet(PlasmaDomain &pd, std::string name);
        
        Grid& grid(int index); 
        Grid& grid(std::string name);
        std::vector<Grid>& allGrids();
        std::string nameFromIndex(int index);
        int indexFromName(std::string name);
        
        int num_variables();
        int num_species();
        bool allGridsInitialized();
        bool allStateGridsInitialized();
        void setOutputFlag(int index, bool flag);
        void setOutputFlag(std::string name, bool flag);
        bool getOutputFlag(int index);
        bool getOutputFlag(std::string name);

        virtual std::vector<std::string> def_var_names() const { 
            assert(false && "Please override def_var_names in derived EquationSet class to define m_var_names");
            return {};
        }

        virtual std::vector<int> state_variables() = 0;

        //Vector quantities (momenta and fields) should be formatted as
        //e.g. { {mom1_x,mom1_y},{mom2_x,mom2_y},etc. }
        //If multiple species are involved, the order in which the indices
        //are returned for each category should match s.t. a given index
        //corresponds to a single species
        //e.g. {density1,density2,etc.} {thermal_energy1,thermal_energy2,etc.}
        //where each number corresponds to a single species
        //If a given quantity is used for multiple species, that index
        //should be returned once for each of these species
        virtual std::vector<int> densities() = 0;
        virtual std::vector<std::vector<int>> momenta() = 0;
        virtual std::vector<int> thermal_energies() = 0;
        virtual std::vector<std::vector<int>> fields() = 0;

        std::vector<Grid> computeTimeDerivatives(double visc_coeff) { return computeTimeDerivatives(m_grids,visc_coeff); }
        void applyTimeDerivatives(const std::vector<Grid> &time_derivatives, double step) { applyTimeDerivatives(m_grids,time_derivatives,step); }

        virtual Grid getDT() = 0;
        virtual std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) = 0;
        virtual void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) = 0;
        virtual void propagateChanges() = 0;
        virtual void computeConstantTerms() = 0;

    protected:
        PlasmaDomain& m_pd;
        std::vector<Grid> m_grids;
        std::vector<bool> m_output_flags;
        const std::vector<std::string> m_var_names;
        std::unordered_map<std::string,int> m_var_indices;
};

#endif