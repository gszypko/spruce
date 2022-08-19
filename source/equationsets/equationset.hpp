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

class PlasmaDomain;

class EquationSet {
    public:
        // Construction
        EquationSet(PlasmaDomain &pd, std::vector<std::string> var_names);
        virtual ~EquationSet() {}
        static std::unique_ptr<EquationSet> instantiateDefault(PlasmaDomain &pd,const std::string &name);
        static std::unique_ptr<EquationSet> instantiateWithConfig(PlasmaDomain &pd, std::ifstream &in_file, const std::string &name, bool active);
        void configureEquationSet(std::ifstream &in_file);
        static const inline std::vector<std::string> m_sets {"ideal_mhd","ideal_mhd_cons","ideal_mhd_2E","ideal_2F"};
        static bool isEquationSetName(const std::string& name);
        virtual std::vector<std::string> config_names() const {return {};};

        // Getters
        Grid& grid(int index);
        Grid& grid(const std::string& name);
        std::vector<Grid> allGrids() const;
        std::vector<std::string> allNames() const;
        std::string nameFromIndex(int index) const;
        int indexFromName(std::string name);
        
        // The number of different variables tracked by the EquationSet
        int num_variables();
        // The number of species tracked by the EquationSEt
        int num_species();
        // Checks if all Grids being tracked have been initialized (size greater than 1)
        bool allGridsInitialized();
        // Checks if all Grids correspond to state variables have been initialized (size greater than 1)
        bool allStateGridsInitialized();
        // Set the flag of the variable corresponding to index, s.t. flag==true means it is written to mhd.out
        void setOutputFlag(int index, bool flag);
        // Set the flag of the variable corresponding to name, s.t. flag==true means it is written to mhd.out
        void setOutputFlag(std::string name, bool flag);
        // Get the flag of the variable corresponding to index; true means it should be written to mhd.out
        bool getOutputFlag(int index);
        // Get the flag of the variable corresponding to index; true means it should be written to mhd.out
        bool getOutputFlag(std::string name);

        // Defines the names corresponding to all of the variables tracked by the EquationSet
        // The ordering of the names here corresponds to the ordering of variables in m_grids
        virtual std::vector<std::string> def_var_names() const { 
            assert(false && "Please override def_var_names in derived EquationSet class to define m_var_names");
            return {};
        }

        // Defines the indices corresponding to all of the state variables for the EquationSet
        // These are the variables that should be written to/read from any .state file
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
        virtual std::vector<int> fields() = 0;

        std::vector<Grid> computeTimeDerivatives(double visc_coeff) { return computeTimeDerivatives(m_grids,visc_coeff); }
        void applyTimeDerivatives(const std::vector<Grid> &time_derivatives, double step) { applyTimeDerivatives(m_grids,time_derivatives,step); }
        void propagateChanges() { propagateChanges(m_grids); }
        void populateVariablesFromState() { populateVariablesFromState(m_grids); }

        // Return a Grid containing a maximum time step for each cell
        virtual Grid getDT() = 0;
        // Return a vector of Grids containing the time derivatives computed for all
        // evolved quantities for which analytic time derivatives exist in the EquationSet. 
        // This vector<Grid> is fed into
        // applyTimeDerivatives() to apply the time evolution to the system.
        // The argument grids should be a vector<Grid> of the same size and dimensions as
        // the member variable m_grids (this allows for evolution of intermediate steps
        // by the PlasmaDomain, without touching the current "real" state of the system)
        virtual std::vector<Grid> computeTimeDerivatives(const std::vector<Grid> &grids, double visc_coeff) = 0;
        // Apply the time evolution described by the result of computeTimeDerivatives()
        // The argument grids should be a vector<Grid> of the same size and dimensions as
        // the member variable m_grids (this allows for evolution of intermediate steps
        // by the PlasmaDomain, without touching the current "real" state of the system)
        virtual void applyTimeDerivatives(std::vector<Grid> &grids, const std::vector<Grid> &time_derivatives, double step) = 0;
        // Make any changes to quantities derived from the evolved quantities to ensure
        // that all quantities are consistent with one another.
        // The argument grids should be a vector<Grid> of the same size and dimensions as
        // the member variable m_grids (this allows for evolution of intermediate steps
        // by the PlasmaDomain, without touching the current "real" state of the system)
        virtual void propagateChanges(std::vector<Grid> &grids) = 0;
        // Populate all tracked quantities, starting from the set of state quantities
        // i.e. quantities read in from the .state file and defined by state_variables()
        // The argument grids should be a vector<Grid> of the same size and dimensions as
        // the member variable m_grids (this allows for evolution of intermediate steps
        // by the PlasmaDomain, without touching the current "real" state of the system)
        virtual void populateVariablesFromState(std::vector<Grid> &grids) = 0;
    protected:
        PlasmaDomain& m_pd;
        std::vector<Grid> m_grids{};
        std::vector<bool> m_output_flags{};
        const std::vector<std::string> m_var_names{};
        std::unordered_map<std::string,int> m_var_indices{};
        virtual void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) = 0;
};

#endif