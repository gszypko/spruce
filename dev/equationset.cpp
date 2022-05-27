#include "equationset.hpp"
#include <iostream>
#include <algorithm>

EquationSet::EquationSet(PlasmaDomain &pd): m_pd(pd) {}

Grid& EquationSet::grid(int index) {
    return m_grids[index];
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
