#include "equationset.hpp"
#include "idealmhd.hpp"
#include "idealmhdcons.hpp"
#include "idealmhd2E.hpp"
#include "ideal2F.hpp"
#include "plasmadomain.hpp"
#include "utils.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <stdexcept>
#include <fstream>


EquationSet::EquationSet(PlasmaDomain &pd, std::vector<std::string> var_names): 
    m_pd(pd),
    m_var_names{var_names}
{
    m_grids = std::vector<Grid>(num_variables(),Grid::Zero(1,1));
    m_output_flags = std::vector<bool>(num_variables(),false);
    for(int i=0; i<m_var_names.size(); i++) m_var_indices[m_var_names[i]] = i;
}

std::unique_ptr<EquationSet> EquationSet::instantiateDefault(PlasmaDomain &pd,const std::string &name)
{
    std::unique_ptr<EquationSet> ptr;
    if (name == "ideal_mhd") ptr = std::unique_ptr<EquationSet>(new IdealMHD(pd));
    else if (name == "ideal_mhd_cons") ptr = std::unique_ptr<EquationSet>(new IdealMHDCons(pd));
    else if (name == "ideal_mhd_2E") ptr = std::unique_ptr<EquationSet>(new IdealMHD2E(pd));
    else if (name == "ideal_2F") ptr = std::unique_ptr<EquationSet>(new Ideal2F(pd));
    else{
        std::cerr << "Equation set name <" << name << "> not recognized." << std::endl;
        assert(false);
    }
    return ptr;
}

void EquationSet::instantiateWithConfig(std::unique_ptr<EquationSet>& eqs,PlasmaDomain &pd, std::ifstream &in_file, const std::string &name, bool active)
{
    if(active){ // instantiate equation set if it is active
        eqs = instantiateDefault(pd,name);
        eqs->configureEquationSet(in_file);    
    }
    else{ // do not instantiate if equation set is not active
        std::string line;
        std::getline(in_file, line);
        clearWhitespace(line);
        assert(line[0] == '{' && "All Equation set activation/deactivation configs must be immediately followed by curly brackets");
        std::getline(in_file, line);
        clearWhitespace(line);
        while(line[0] != '}'){
            std::getline(in_file, line);
            clearWhitespace(line);
        }
    }
}

bool EquationSet::isEquationSetName(const std::string& name)
{
    auto it = std::find(m_sets.begin(),m_sets.end(),name);
    return it != m_sets.end();
}

void EquationSet::configureEquationSet(std::ifstream &in_file)
{
    std::vector<std::string> lhs_strings, rhs_strings;
    std::string line;
    std::getline(in_file, line);
    assert(line[0] == '{' && "All equation set activation configs must be immediately followed by curly brackets (on their own lines) to enclose Module configs");
    if (line[0] == '{'){
        std::getline(in_file, line);
        while(line[0] != '}'){
            clearWhitespace(line);
            if (line.empty() || line[0] == '#') continue; // skip comment and empty lines
            std::istringstream ss_line(line);
            std::string lhs, rhs;
            std::getline(ss_line,lhs,'=');
            std::getline(ss_line,rhs,'#');
            lhs_strings.push_back(lhs);
            rhs_strings.push_back(rhs);
            std::getline(in_file, line);
        }
    }
    parseEquationSetConfigs(lhs_strings,rhs_strings);
}

void EquationSet::setupEquationSet()
{
    assert(allStateGridsInitialized() && "All variables specified as state variables for the current EquationSet must be specified in the .state file");
    populateVariablesFromState();
    assert(allGridsInitialized() && "All variables must be initialized by <populateVariablesFromState>.");
    m_grids_dt = std::vector<Grid>(num_variables(),Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    setupEquationSetDerived();
    indexFromName("dt"); // ensure that dt is a grid within the equation set, this function throws error if not found
}

void EquationSet::populateVariablesFromState(std::vector<Grid> &grids){
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    recomputeEvolvedVarsFromStateVars(grids);
    m_pd.updateGhostZones();
    recomputeDerivedVarsFromEvolvedVars(grids);
    recomputeDT();
}

// any steps required for setup after instantiation of module but before iteration of grids
// typical use is to instantiate sizes of internal grids and compute any grids that are constant in time
void EquationSet::setupEquationSetDerived()
{
    return;
}

Grid& EquationSet::grid(int index) {
    assert(index >= 0 && index < num_variables() && "Grid index must be within range of m_grids");
    return m_grids[index];
}

Grid& EquationSet::grid(const std::string& name) {
    int index = indexFromName(name);
    return grid(index);
}

bool EquationSet::allGridsInitialized() {
    for(const Grid& g : m_grids){
        if(g.size() == 1) return false;
    }
    return true;
}

bool EquationSet::allStateGridsInitialized() {
    for(int i : state_variables()){
        if(grid(i).size() == 1) return false;
    }
    return true;
}

std::vector<Grid> EquationSet::allGrids() const {
    return m_grids;
}

std::vector<std::string> EquationSet::allNames() const {
    return m_var_names;
}

std::string EquationSet::nameFromIndex(int index) const {
    return m_var_names[index];
}

int EquationSet::num_variables() {
    return m_var_names.size();
}

int EquationSet::num_species() {
    assert(densities().size() == momenta().size() && momenta().size() == thermal_energies().size() && \
            "Each species must have an entry each in densities, momenta, thermal_energies, and fields");
    return densities().size();
}

void EquationSet::setOutputFlag(int index, bool flag){
    m_output_flags[index] = flag;
}

void EquationSet::setOutputFlag(std::string name, bool flag){
    m_output_flags[indexFromName(name)] = flag;
}

bool EquationSet::getOutputFlag(int index){
    return m_output_flags[index];
}

bool EquationSet::getOutputFlag(std::string name){
    return m_output_flags[indexFromName(name)];
}

int EquationSet::indexFromName(std::string name){
    int ind{};
    try { ind = m_var_indices.at(name); }
    catch (const std::out_of_range& e) {
        std::cerr << "Variable name <" << name << "> not recognized\n";
        assert(false);
        return 0;
    }
    return ind;
}

Grid EquationSet::getDT() {return m_grids[indexFromName("dt")];}

void EquationSet::applyTimeDerivativesBase(std::vector<Grid> &grids, double step)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    for (int i : evolved_variables()) grids[i] += step*m_grids_dt[i];
    propagateChanges(grids);
}

void EquationSet::applyTimeDerivativesBase(std::vector<Grid> &grids,const std::vector<Grid> &grids_dt, double step)
{
    assert(grids.size() == m_grids.size() && "This function designed to operate on full system vector<Grid>");
    for (int i : evolved_variables()) grids[i] += step*grids_dt[i];
    propagateChanges(grids);
}