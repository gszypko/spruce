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
        // ***************** STATIC MEMBERS

        static const inline std::vector<std::string> m_sets {"ideal_mhd","ideal_mhd_cons","ideal_mhd_2E","ideal_2F"};
        static bool isEquationSetName(const std::string& name);

        // ***************** CONSTRUCTION - to be called while reading .config file

        EquationSet(PlasmaDomain &pd, std::vector<std::string> var_names);
        virtual ~EquationSet() {}
        static std::unique_ptr<EquationSet> instantiateDefault(PlasmaDomain &pd,const std::string &name);
        static void instantiateWithConfig(std::unique_ptr<EquationSet>& eqs,PlasmaDomain &pd, std::ifstream &in_file, const std::string &name, bool active);
        void configureEquationSet(std::ifstream &in_file);
        
        // ***************** SETUP FUNCTIONS - to be called during PlasmaDomain construction after the .state file has been read

        void setupEquationSet();
        virtual void setupEquationSetDerived();
        virtual std::vector<std::string> config_names() const {return {};};

        // ***************** VARIABLE DEFINITIONS 

        // Defines the names corresponding to m_grids
        virtual std::vector<std::string> def_var_names() const { 
            assert(false && "Please override def_var_names in derived EquationSet class to define m_var_names");
            return {};
        }
        // Defines the indices corresponding to the state variables
        virtual std::vector<int> state_variables() const = 0;
        // Defines the indices corresponding to all of the evolved variables
        virtual std::vector<int> evolved_variables() const = 0;

        // ***************** VARIABLE INDICES FOR EVOLVED VARIABLES

        // indices corresponding to the evolved density of each plasma species (i.e., {i_rho,e_rho} for a 2F model)
        virtual std::vector<std::string> species() const = 0;
        // indices corresponding to the evolved density of each plasma species (i.e., {i_rho,e_rho} for a 2F model)
        virtual std::vector<int> densities() const = 0;
        // indices corresponding to the evolved density of each plasma species (i.e., {i_rho,e_rho} for a 2F model)
        virtual std::vector<int> number_densities() const = 0;
        // vector of indices to the evolved momenta for each plasma species (i.e., {{i_mom_x,i_mom_y},{e_mom_x,e_mom_y}} for a 2F model)
        virtual std::vector<std::vector<int>> momenta() const = 0;
        // vector of indices to the velocity for each plasma species (i.e., {{i_v_x,i_v_y},{e_v_x,e_v_y}} for a 2F model)
        virtual std::vector<std::vector<int>> velocities() const = 0;
        // indices corresponding to the evolved thermal energy of each plasma species (i.e., {i_thermal_energy,e_thermal_energy} for a 2F model)
        virtual std::vector<int> thermal_energies() const = 0;
        // indices corresponding to the evolved pressures of each plasma species (i.e., {i_thermal_energy,e_thermal_energy} for a 2F model)
        virtual std::vector<int> pressures() const = 0;
        // indices corresponding to the evolved thermal energy of each plasma species (i.e., {i_thermal_energy,e_thermal_energy} for a 2F model)
        virtual std::vector<int> temperatures() const = 0;
        // indices corresponding to the evolved electromagnetic fields
        virtual std::vector<int> fields() const = 0;
        // indices corresponding to the dt grid for each species
        virtual std::vector<int> timescale() const = 0;

        // ***************** NON-CONST SETTERS

        // Return grid corresponding to index
        Grid& grid(int index);
        // Return grid corresponding to grid name
        Grid& grid(const std::string& name);
        // Set the flag of the variable corresponding to index, s.t. flag==true means it is written to mhd.out
        void setOutputFlag(int index, bool flag);
        // Set the flag of the variable corresponding to name, s.t. flag==true means it is written to mhd.out
        void setOutputFlag(std::string name, bool flag);

        // ***************** CONST GETTERS

        // Returns m_grids
        std::vector<Grid> allGrids() const;
        // Returns m_grids_dt
        std::vector<Grid> allGridsDT() const;
        // Returns m_var_names
        std::vector<std::string> allNames() const;
        // Returns grid name corresponding to index
        std::string index2name(int index) const;
        // Returns grid index corresonding to name
        int name2index(std::string name) const;
        // Returns the evolved variable index for given grid name
        int name2evolvedindex(std::string name) const;
        // Checks whether <name> is a valid grid name
        bool is_var(std::string name) const;
        // Checks whether <name> is a valid state variable
        bool is_state_var(std::string name) const;
        // Checks whether <name> is a valid evolved variable
        bool is_evolved_var(std::string name) const;
        // Return m_grids[dt]
        Grid getDT() const;
        // The number of different variables tracked by the EquationSet
        int num_variables() const;
        // The number of species tracked by the EquationSEt
        int num_species() const;
        // Checks if all Grids being tracked have been initialized (size greater than 1)
        bool allGridsInitialized() const;
        // Checks if all Grids correspond to state variables have been initialized (size greater than 1)
        bool allStateGridsInitialized() const;
        // Get the flag of the variable corresponding to index; true means it should be written to mhd.out
        bool getOutputFlag(int index) const;
        // Get the flag of the variable corresponding to index; true means it should be written to mhd.out
        bool getOutputFlag(std::string name) const;

        
        // ***************** FUNCTIONS FOR EVOLVING VARIABLES
        
        // Return time derivatives for evolved variables - calls derived class and module compute time derivative functions
        std::vector<Grid> computeTimeDerivatives(const std::vector<Grid>& grids) const;
        // Overload for computing time derivatives for internal m_grids
        std::vector<Grid> computeTimeDerivatives() const { return computeTimeDerivatives(m_grids); };
        // Apply the time derivatives to the evolved variables and update derived grids
        void applyTimeDerivatives(std::vector<Grid>& grids, const std::vector<Grid>& grids_dt, double step) const;
        // Overload for applying time derivatives to m_grids
        void applyTimeDerivatives(const std::vector<Grid>& grids_dt,double step) { applyTimeDerivatives(m_grids,grids_dt,step); };
        // Updated derived variables - to be called after time derivatives are applied to evolved variables
        void propagateChanges(std::vector<Grid>& grids) const;
        // Overload for propagating changes to m_grids
        void propagateChanges() { propagateChanges(m_grids); };
        // enforce minimums to grids
        virtual void enforceMinimums(std::vector<Grid>& grids) const = 0;
        // enforce minimums to grids
        void enforceMinimums() {enforceMinimums(m_grids); };
        
    protected:
        // reference to the top-level simulation class
        PlasmaDomain& m_pd; // NOTE: should this be const?

        // full vector of equation set grids, their names, corresponding vector indices
        std::vector<Grid> m_grids{};
        const std::vector<std::string> m_var_names{};
        std::unordered_map<std::string,int> m_var_indices{};
        std::vector<bool> m_output_flags{}; // indicates whether a given variable should be written to the .out file

        // ***************** PRIVATE MEMBER FUNCTIONS     
        // Parse configs to initialize class members
        virtual void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) = 0;
        // Derived class function for returning time derivatives of evolved variables
        virtual std::vector<Grid> computeTimeDerivativesDerived(const std::vector<Grid>& grids) const = 0;
        // Populate evolved/derived variables from state variables
        void populateVariablesFromState(std::vector<Grid>& grids) const;
        // Overload to apply populate variables for internal m_grids
        void populateVariablesFromState() { populateVariablesFromState(m_grids); };
        // Compute evolved variables from the state variables        
        virtual void recomputeEvolvedVarsFromStateVars(std::vector<Grid>& grids) const = 0;
        // Overload to compute evolved variables using internal m_grids
        void recomputeEvolvedVarsFromStateVars() { recomputeEvolvedVarsFromStateVars(m_grids); };
        // Recompute derived variables from the evolved variables
        virtual void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid>& grids) const = 0;
        // Overload to recompute derived variables for internal m_grids
        void recomputeDerivedVarsFromEvolvedVars() { recomputeDerivedVarsFromEvolvedVars(m_grids); };
        // Recompute grids[dt]
        virtual void recomputeDT(std::vector<Grid>& grids) const = 0;
        // Apply recomputeDT to internal m_grids
        void recomputeDT() { recomputeDT(m_grids); };
};

#endif