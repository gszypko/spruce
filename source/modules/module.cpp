//module.cpp

#include "module.hpp"
#include <vector>

Module::Module(PlasmaDomain &pd, bool propagate_changes): m_pd(pd), m_propagate_changes(propagate_changes) {}

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
            if(line.empty() || line[0] == '#') continue; //skip comment and empty lines
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

std::string Module::commandLineMessage() const
{
    return "";
}

void Module::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const
{
    return;
}
