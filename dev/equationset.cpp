#include "equationset.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>

EquationSet::EquationSet(PlasmaDomain &pd, std::vector<std::string> var_names): 
    m_pd(pd),
    m_var_names{var_names}
{
    m_grids = std::vector<Grid>(num_variables(),Grid::Zero(1,1));
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
    auto it = std::find(m_var_names.begin(),m_var_names.end(),name);
    if(it == m_var_names.end()){
        std::cerr << "Variable name " << name << " not recognized\n";
        abort();
    }
    auto index = std::distance(m_var_names.begin(),it);
    return grid(index);
}

std::string EquationSet::gridName(int index) {
    return m_var_names[index];
}

int EquationSet::num_variables() {
    return m_var_names.size();
}
