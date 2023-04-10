//module.cpp

#include "module.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <cassert>

Module::Module(PlasmaDomain &pd): m_pd(pd) {}

void Module::configureModule(std::ifstream &in_file)
{
    std::vector<std::string> lhs_strings, rhs_strings;
    std::string line;
    std::getline(in_file, line);
    assert(line[0] == '{' && "All Module activation configs must be immediately followed by curly brackets (on their own lines) to enclose Module configs");
    if(line[0] == '{'){
        std::getline(in_file, line);
        while(line[0] != '}'){
            clearWhitespace(line);
            if(line.empty() || line[0] == '#'){ //skip comment and empty lines
                std::getline(in_file, line);
                continue; 
            }
            std::istringstream ss_line(line);
            std::string lhs, rhs;
            std::getline(ss_line,lhs,'=');
            std::getline(ss_line,rhs,'#');
            lhs_strings.push_back(lhs);
            rhs_strings.push_back(rhs);
            std::getline(in_file, line);
        }
    }
    parseModuleConfigs(lhs_strings,rhs_strings);
}

// any steps required for setup after instantiation of module but before iteration of grids
// typical use is to instantiate sizes of internal grids and compute any grids that are constant in time
void Module::setupModule()
{
    return;
}

void Module::iterateModule(double dt)
{
    return;
}

void Module::preIterateModule(double dt)
{
    return;
}

void Module::postIterateModule(double dt)
{
    return;
}

void Module::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt)
{
    return;
}

void Module::preRecomputeDerivedModule(std::vector<Grid>& grids) const
{
    return;
}

std::string Module::commandLineMessage() const
{
    return "";
}

void Module::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    return;
}
