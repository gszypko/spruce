//modulehandler.hpp
//Defines ModuleHandler class to manage Modules of behavior
//on top of base MHD evolution in a PlasmaDomain

#ifndef MODULEHANDLER_HPP
#define MODULEHANDLER_HPP

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include "grid.hpp"
#include "module.hpp"

class PlasmaDomain;

class ModuleHandler {
    public:
        ModuleHandler(PlasmaDomain &pd);
        void setupModules();
        void iterateModules(double dt);
        void preIterateModules(double dt);
        void postIterateModules(double dt);
        void iterateComputeTimeDerivativesModules(const std::vector<Grid> &grids,std::vector<Grid>& grids_dt);
        void iteratePreRecomputeDerivedModule(std::vector<Grid>& grids) const;
        void instantiateModule(const std::string &module_name, std::ifstream &in_file, bool module_active = true);
        bool isModuleName(std::string name) const; //Check if string is valid module name defined within this class
        std::vector<std::string> getCommandLineMessages() const; //Return short messages from Modules to print to command line
        void getFileOutputData(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const; //Get variable names and Grids to print to output file

    private:
        PlasmaDomain &m_pd;
        std::vector<std::unique_ptr<Module>> m_modules;
        static inline std::vector<std::string> m_module_names = {
            "radiative_losses", "thermal_conduction", "ambient_heating", "anomalous_resistivity", "localized_heating", "field_heating", "sg_filtering",
            "coulomb_explosion","eic_thermalization", "artificial_viscosity", "global_temperature"
        };
};

#endif