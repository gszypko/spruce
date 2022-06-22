#ifndef PLASMA_SETTINGS_HPP
#define PLASMA_SETTINGS_HPP

#include "settings.hpp"
#include <unordered_map>

class PlasmaSettings : public Settings
{
public: 
    // *** Construction
    PlasmaSettings(fs::path settings_path,std::string unit = "cgs");
    void update_possible_units();
    std::string dependencies(std::string characteristic);
private:
    // *** Members
    const str_vec m_plasma_characteristics {"w_pi","w_pe","l_deb","sig","tau_exp"};
    std::unordered_map<std::string,int> m_ind;
};

#endif