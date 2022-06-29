#include "ConfigUI.hpp"
#include <fstream>
#include <cassert>
#include "eic_thermalization.hpp"
#include "coulomb_explosion.hpp"
#include <sstream>

ConfigUI::ConfigUI(PlasmaDomain& pd, fs::path config_path): 
    m_config_path{config_path}, m_pd{pd}
{
    instantiate_modules();
    define_config_names();
    read_config();
}

void ConfigUI::read_config()
{
    assert(m_config_path.extension().string()==".config");
    assert(fs::exists(m_config_path) && !fs::is_empty(m_config_path) && "File must exist and not be empty.");
    m_config_data.reserve(100);
    std::ifstream file_stream(m_config_path);
    while (file_stream.good()){
        std::string curr_line;
        std::getline(file_stream,curr_line);
        m_config_data.push_back(curr_line);
        if (m_config_data.size() == m_config_data.capacity()) 
            m_config_data.reserve(2*m_config_data.capacity());
    }
    file_stream.close();
    m_config_data.shrink_to_fit();
}

void ConfigUI::define_config_names()
{
    for (auto i : m_additional_configs)
        m_config_names.push_back(i);
    for (const auto& module : m_modules)
        for (auto i : module->config_names())
            m_config_names.push_back(i);
}

std::string ConfigUI::get_eqnset_name() const
{
    for (auto line : m_config_data){
        std::string lhs,rhs;
        std::istringstream ss(line);
        std::getline(ss,lhs,'=');
        std::getline(ss,rhs);
        std::string lhs_cleared = lhs; lhs_cleared.erase(std::remove(lhs_cleared.begin(),lhs_cleared.end(),' '),lhs_cleared.end());
        if (lhs_cleared == "equation_set") return rhs;
    }
    assert(false && "Equation set not found in .config file.");
    return "";
}

void ConfigUI::instantiate_modules()
{
    for (auto name : module_names){
        if (name == "eic_thermalization") m_modules.push_back(std::unique_ptr<Module>(new EICThermalization(m_pd)));
        else if (name == "coulomb_explosion") m_modules.push_back(std::unique_ptr<Module>(new CoulombExplosion(m_pd)));
        else assert(false && "Module name not recognized.");
    }
}

void ConfigUI::update_config(std::string name,std::string val)
{
    // ensure that input is a config name
    auto it = std::find(m_config_names.begin(),m_config_names.end(),name);
    assert(it!=m_config_names.end());
    // find the line in config data that starts with this value
    bool config_found{false};
    for (int i=0; i<m_config_data.size(); i++){
        std::string curr_line = m_config_data[i];
        std::string lhs,rhs;
        std::istringstream ss(curr_line);
        std::getline(ss,lhs,'=');
        std::getline(ss,rhs);
        std::string lhs_cleared = lhs; lhs_cleared.erase(std::remove(lhs_cleared.begin(),lhs_cleared.end(),' '),lhs_cleared.end());
        if (lhs_cleared.compare(name)==0){
             m_config_data[i] = lhs + " = " + val;
             config_found = true;
        }
    }
    if (!config_found) m_config_data.push_back(name + " = " + val);
    
}

bool ConfigUI::is_config(std::string name) const
{
    auto it = std::find(m_config_names.begin(),m_config_names.end(),name);
    return it != m_config_names.end();
}

void ConfigUI::write_config_file(const fs::path& directory) const
{
    if (!fs::exists(directory)) fs::create_directories(directory);
    fs::path filename = "ucnp.config";
    std::ofstream outfile(directory/filename);
    for (const auto& curr_line : m_config_data){
        outfile << curr_line << std::endl;
    }
    outfile.close();
}