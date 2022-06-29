#ifndef CONFIG_UI_HPP
#define CONFIG_UI_HPP

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <module.hpp>

class PlasmaDomain;

namespace fs = std::filesystem;

class ConfigUI
{
public:   
    // *** Manual Controls
    const std::vector<std::string> module_names {"eic_thermalization","coulomb_explosion"};
    const std::vector<std::string> m_additional_configs {"duration","time_output_interval"};
    // *** Construction
    ConfigUI(PlasmaDomain& pd, fs::path config_path);
    ConfigUI()=default;
    void read_config();
    void instantiate_modules();
    void define_config_names();
    // *** Usage
    void update_config(std::string name,std::string val);
    void write_config_file(const fs::path& directory) const;
    bool is_config(std::string name) const;
    std::string get_eqnset_name() const;
private:
    // *** Members
    fs::path m_config_path{};
    std::vector<std::string> m_config_data{};
    std::vector<std::string> m_config_names{};
    std::vector<std::unique_ptr<Module>> m_modules{};
    PlasmaDomain& m_pd;
};

#endif