#ifndef UCNP_SETTINGS_HPP
#define UCNP_SETTINGS_HPP

#include "settings.hpp"
#include <unordered_map>

class UCNP : public Settings
{
public: 
    // *** Construction
    UCNP(fs::path settings_path,std::string unit = "cgs");
    void update_possible_units() override;
    str_vec get_dependencies(std::string plasma_characteristic) const;
    void check_for_dependencies() const;
    double get_characteristic(std::string str) const;
    void update_variables();
private:
    // *** Members
    const str_vec m_plasma_characteristics {"w_pi","w_pe","l_deb","sig","tau_exp","a"};
    enum Vars {w_pi,w_pe,l_deb,sig,tau_exp,a};
};

#endif