#ifndef STATE_HANDLER_HPP
#define STATE_HANDLER_HPP

#include <iostream>
#include <filesystem>
#include <vector>
#include <string>
#include <unordered_map>
#include "plasmadomain.hpp"
#include "equationset.hpp"
#include "grid.hpp"
#include "settings.hpp"

namespace fs = std::filesystem;

class StateHandler
{
public:
    // *** Construction
    StateHandler(const std::string& eqs_name);
    // *** Usage
    int gridname2ind(const std::string& name) const;
    int varname2ind(const std::string& name) const;
    void setvar(const std::string& name, double val);
    void setgrid(const std::string& name, Grid grid);
    double getvar(const std::string& name) const;
    Grid getgrid(const std::string& name) const;
    void all_initialized() const;
    void write_state_file(const fs::path& directory) const;
    // *** Equation Set Initializers
    void setup(const std::unique_ptr<Settings>& pms);
    void setup_idealmhd(const std::unique_ptr<Settings>& pms);
    void setup_idealmhd2e(const std::unique_ptr<Settings>& pms);
private:
    // *** Members - Equation Set
    PlasmaDomain m_pd{};
    std::string m_eqs_name;
    std::unique_ptr<EquationSet> m_eqs;
    // *** Members - Preamble Variables
    std::vector<std::string> m_varnames {"xdim","ydim","ion_mass","adiabatic_index","time"};
    std::vector<double> m_vars {std::vector(m_varnames.size(),0.0)};
    std::vector<bool> m_vars_initialized {std::vector(m_varnames.size(),false)};
    // *** Members - Grids
    std::vector<std::string> m_gridnames;
    std::vector<Grid> m_grids;
    std::vector<bool> m_grids_initialized;
};

#endif