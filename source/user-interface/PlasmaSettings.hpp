#ifndef PLASMA_SETTINGS_HPP
#define PLASMA_SETTINGS_HPP

#include "ui_utility.hpp"

#include <iostream> // std::cerr, std::endl, std::cout
#include <vector> // std::vector
#include <string> // std::string, std::getline
#include <filesystem> // std::filesystem
namespace fs = std::filesystem;
#include <fstream> // std::ofstream, std::ifstream
#include <sstream> // std::istringstream

class PlasmaSettings
{
public: 
    // constructors

    PlasmaSettings(fs::path settings_path, int array);
    // PlasmaSettings& operator=(const PlasmaSettings& other);
    
    // member functions

    void load_settings(const fs::path& settings_path);
    void check_units(void) const;
    double getvar(const std::string& name) const; 
    std::string getopt(const std::string& name) const;
    void write_array_params(const fs::path& path) const;
private:
    // member variables corresponding to specified m_array
    const int m_array;
    std::vector<std::string> m_names;
    std::vector<std::string> m_units;
    std::vector<std::string> m_val; // m_val = m_unique_val[array]

    // matrix of unique simulation parameters
    std::vector<std::vector<std::string>> m_mat; // m_mat[i] is a vector of values corresponding to m_names in terms of m_units
};

#endif