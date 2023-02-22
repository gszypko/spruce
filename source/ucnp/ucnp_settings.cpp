#include "ucnp_settings.hpp"
#include "ucnputils.hpp"
#include <cmath>
#include "constants.hpp"

UCNP::UCNP(fs::path settings_path,std::string unit): Settings(unit)
{
    load_settings(settings_path);
    if (m_runs_found) process_runs();
    m_possible_units = get_units_from_vars();
    update_possible_units();
    check_units();
    choose_array();
    check_for_dependencies();
    update_variables();
    choose_array();
}

// append the plasma characteristics to the list of possible units
void UCNP::update_possible_units()
{
    for (auto i : m_plasma_characteristics) m_possible_units.push_back(i);
}

// return the dependencies for each plasma characteristic
Settings::str_vec UCNP::get_dependencies(std::string str) const
{
    auto it = std::find(m_plasma_characteristics.begin(),m_plasma_characteristics.end(),str);
    assert(it!=m_plasma_characteristics.end());
    int loc = std::distance(m_plasma_characteristics.begin(),it);
    switch (loc)
    {   case w_pi: return {"n","m_i"};
        case w_pe: return {"n"};
        case l_deb: return {"n","Te"};
        case sig: return {"sig_x","sig_y"};
        case tau: return {"m_i","Te","Ti","sig_x","sig_y"};
        case tau_x: return {"m_i","Te","Ti","sig_x"};
        case tau_y: return {"m_i","Te","Ti","sig_y"};
        case a: return {"n"};
        default: return {""};
    }
}

// ensure that all dependencies are found within .settings file
void UCNP::check_for_dependencies() const
{
    for (int i=0; i<m_plasma_characteristics.size(); i++){
        str_vec dependencies = get_dependencies(m_plasma_characteristics[i]);
        for (int j=0; j<dependencies.size(); j++){
            if (!is_name(dependencies[j])){
                std::cerr << "Dependency <" << dependencies[j] << "> is not a possible unit." << std::endl;
                assert(false);
            }
        }
    }
}

double UCNP::get_characteristic(std::string str) const
{
    double result;
    auto it = std::find(m_plasma_characteristics.begin(),m_plasma_characteristics.end(),str);
    assert(it!=m_plasma_characteristics.end());
    int loc = std::distance(m_plasma_characteristics.begin(),it);
    switch (loc)
    {   case w_pi: return phys::plasma_freq(getval("n"),getval("m_i"));
        case w_pe: return phys::plasma_freq(getval("n"),E);
        case l_deb: return phys::debye_length(getval("n"),getval("Te"));
        case sig: return pow(getval("sig_x")*getval("sig_y"),1./2.);
        case tau: return phys::tau_exp(get_characteristic("sig"),getval("m_i"),getval("Te")+getval("Ti"));
        case tau_x: return phys::tau_exp(getval("sig_x"),getval("m_i"),getval("Te")+getval("Ti"));
        case tau_y: return phys::tau_exp(getval("sig_y"),getval("m_i"),getval("Te")+getval("Ti"));
        case a: return phys::wigner_seitz_radius(getval("n"));
        default: return -1.;
    }
}

void UCNP::update_variables()
{
    // for each plasma characteristic, append the characteristic name and unit
    // update each row of m_unique with the values of each characteristic
    for (int i=0; i<m_plasma_characteristics.size(); i++){
        m_names.push_back(m_plasma_characteristics[i]);
        m_units.push_back(m_unit_str);
        for (int j=0; j<m_unique.size(); j++){
            m_unique[j].push_back(num2str(get_characteristic(m_plasma_characteristics[i])));
        }
    }
}