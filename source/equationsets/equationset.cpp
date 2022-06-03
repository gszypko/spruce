#include "equationset.hpp"
#include "idealmhd.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <stdexcept>

EquationSet::EquationSet(PlasmaDomain &pd, std::vector<std::string> var_names): 
    m_pd(pd),
    m_var_names{var_names}
{
    m_grids = std::vector<Grid>(num_variables(),Grid::Zero(1,1));
    m_output_flags = std::vector<bool>(num_variables(),false);
    for(int i=0; i<m_var_names.size(); i++) m_var_indices[m_var_names[i]] = i;
}

std::unique_ptr<EquationSet> EquationSet::spawnEquationSet(PlasmaDomain &pd, std::string name){
    if(name == "ideal_mhd") return std::unique_ptr<EquationSet>(new IdealMHD(pd));
    else{
        std::cerr << "EquationSet name " << name << " not recognized\n";
        abort();
    }
}

Grid& EquationSet::grid(int index) {
    assert(index >= 0 && index < num_variables() && "Grid index must be within range of m_grids");
    return m_grids[index];
}

bool EquationSet::allGridsInitialized() {
    for(Grid g : m_grids){
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

Grid& EquationSet::grid(std::string name) {
    int index = indexFromName(name);
    return grid(index);
}

std::vector<Grid>& EquationSet::allGrids() {
    return m_grids;
}

std::string EquationSet::nameFromIndex(int index) {
    return m_var_names[index];
}

int EquationSet::num_variables() {
    return m_var_names.size();
}

int EquationSet::num_species() {
    assert(densities().size() == momenta().size() && momenta().size() == thermal_energies().size() \
            && thermal_energies().size() == fields().size() && \
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
    try { return m_var_indices.at(name); }
    catch (const std::out_of_range& e) {
        std::cerr << "Variable name " << name << " not recognized\n";
        assert(false && "Variable name was not recognized");
    }
    // auto it = std::find(m_var_names.begin(),m_var_names.end(),name);
    // if(it == m_var_names.end()){
    //     std::cerr << "Variable name " << name << " not recognized\n";
    //     abort();
    // }
    // auto index = std::distance(m_var_names.begin(),it);
    // return index;
}

