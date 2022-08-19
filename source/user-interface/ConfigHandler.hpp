#ifndef CONFIG_UI_HPP
#define CONFIG_UI_HPP

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <module.hpp>
#include "plasmadomain.hpp"

namespace fs = std::filesystem;

class ConfigHandler
{
public:   
    // *** Construction
    ConfigHandler(fs::path config_path);
    ConfigHandler()=default;
    void read_config();
    void instantiate_modules();
    void instantiate_equation_set();
    void define_config_names();
    // *** Usage
    void update_config(std::string name,std::string val);
    void write_config_file(const fs::path& directory) const;
    bool is_config(std::string name) const;
    void parse_config_line(const std::string& line,std::string& lhs,std::string& rhs) const;
    std::string eqs_set_name() {return m_eqn_set_name;}
private:
    // *** Members
    const std::vector<std::string> module_names {"eic_thermalization","coulomb_explosion"};
    PlasmaDomain m_pd; // defualt-initialized PlasmaDomain class, used to get list of config names
    const fs::path m_config_path{}; // relative or full path to .config file
    std::vector<std::string> m_config_data{}; // each element contains one un-modified line of the .config file
    std::vector<std::string> m_config_names{}; // full list of possible config names (module, PlasmaDomain)
    std::vector<std::unique_ptr<Module>> m_modules{}; // vector of modules, need to be instantiated to get config names
    std::string m_eqn_set_name{}; // name for equation set, obtained from .config file
    std::unique_ptr<EquationSet> m_eqs; // container for equation set, used to get list of configs
};

#endif