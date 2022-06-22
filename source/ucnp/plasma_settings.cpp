#include "plasma_settings.hpp"

PlasmaSettings::PlasmaSettings(fs::path settings_path,std::string unit):
    Settings(settings_path,unit)
{
    define_possible_units();
}

void PlasmaSettings::update_possible_units()
{
    Settings::define_possible_units();
    // for (const auto& name : plasma_characteristics)
}

std::string PlasmaSettings::dependencies(std::string characteristic)
{
    // switch (m_ind.at(characteristic))
    return ";";
}