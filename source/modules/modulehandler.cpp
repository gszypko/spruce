//modulehandler.cpp

#include "utils.hpp"
#include "modulehandler.hpp"
#include "module.hpp"
#include "ambientheating.hpp"
#include "localizedheating.hpp"
#include "fieldheating.hpp"
#include "thermalconduction.hpp"
#include "radiativelosses.hpp"
#include "anomalousresistivity.hpp"
#include "momentuminjection.hpp"
#include "sgfilter.hpp"
#include "coulomb_explosion.hpp"
#include "eic_thermalization.hpp"
#include "viscosity.hpp"
#include "global_temperature.hpp"
#include "tracerparticles.hpp"
#include "massinjection.hpp"
#include "ambientheatingsink.hpp"
#include "physicalviscosity.hpp"
#include "divcleaning.hpp"

ModuleHandler::ModuleHandler(PlasmaDomain &pd): m_pd(pd) {}

void ModuleHandler::setupModules()
{
    for(int i=0; i<m_modules.size(); i++){
        m_modules[i]->setupModule();
    }
}

void ModuleHandler::iterateModules(double dt)
{
    for(int i=0; i<m_modules.size(); i++){
        m_modules[i]->iterateModule(dt);
    }
}

void ModuleHandler::preIterateModules(double dt)
{
    for(int i=0; i<m_modules.size(); i++){
        m_modules[i]->preIterateModule(dt);
    }
}

void ModuleHandler::postIterateModules(double dt)
{
    for(int i=0; i<m_modules.size(); i++){
        m_modules[i]->postIterateModule(dt);
    }
}

void ModuleHandler::iterateComputeTimeDerivativesModules(const std::vector<Grid> &grids,std::vector<Grid>& grids_dt)
{
    for (int i=0; i<m_modules.size(); i++)
        m_modules[i]->computeTimeDerivativesModule(grids,grids_dt);
}

void ModuleHandler::iteratePreRecomputeDerivedModule(std::vector<Grid>& grids) const
{
    for (int i=0; i<m_modules.size(); i++)
        m_modules[i]->preRecomputeDerivedModule(grids);
}

bool ModuleHandler::isModuleName(std::string name) const
{
    auto it = std::find(m_module_names.begin(),m_module_names.end(),name);
    return (it != m_module_names.end());
}

//Note: assumes that module_name has already been determined
//to be a valid Module name (using isModuleName())
//Setting module_active to false will NOT instantiate the Module,
//instead fast-forwarding in_file to the end of the Module configs
void ModuleHandler::instantiateModule(const std::string &module_name, std::ifstream &in_file, bool module_active)
{
    if(!module_active){
        std::string line;
        std::getline(in_file, line);
        clearWhitespace(line);
        assert(line[0] == '{' && "All Modules activation/deactivation configs must be immediately followed by curly brackets");
        std::getline(in_file, line);
        clearWhitespace(line);
        while(line[0] != '}'){
            std::getline(in_file, line);
            clearWhitespace(line);
        }
        return;
    }
    if(module_name == "ambient_heating") m_modules.push_back(std::unique_ptr<Module>(new AmbientHeating(m_pd)));
    else if(module_name == "radiative_losses") m_modules.push_back(std::unique_ptr<Module>(new RadiativeLosses(m_pd)));
    else if(module_name == "thermal_conduction") m_modules.push_back(std::unique_ptr<Module>(new ThermalConduction(m_pd)));
    else if(module_name == "localized_heating") m_modules.push_back(std::unique_ptr<Module>(new LocalizedHeating(m_pd)));
    else if(module_name == "field_heating") m_modules.push_back(std::unique_ptr<Module>(new FieldHeating(m_pd)));
    else if(module_name == "sg_filtering") m_modules.push_back(std::unique_ptr<Module>(new SGFilter(m_pd)));
    else if(module_name == "coulomb_explosion") m_modules.push_back(std::unique_ptr<Module>(new CoulombExplosion(m_pd)));
    else if(module_name == "eic_thermalization") m_modules.push_back(std::unique_ptr<Module>(new EICThermalization(m_pd)));
    else if(module_name == "artificial_viscosity") m_modules.push_back(std::unique_ptr<Module>(new Viscosity(m_pd)));
    else if(module_name == "global_temperature") m_modules.push_back(std::unique_ptr<Module>(new GlobalTemperature(m_pd)));
    else if(module_name == "anomalous_resistivity") m_modules.push_back(std::unique_ptr<Module>(new AnomalousResistivity(m_pd)));
    else if(module_name == "momentum_injection") m_modules.push_back(std::unique_ptr<Module>(new MomentumInjection(m_pd)));
    else if(module_name == "tracer_particles") m_modules.push_back(std::unique_ptr<Module>(new TracerParticles(m_pd)));
    else if(module_name == "mass_injection") m_modules.push_back(std::unique_ptr<Module>(new MassInjection(m_pd)));
    else if(module_name == "ambient_heating_sink") m_modules.push_back(std::unique_ptr<Module>(new AmbientHeatingSink(m_pd)));
    else if(module_name == "physical_viscosity") m_modules.push_back(std::unique_ptr<Module>(new PhysicalViscosity(m_pd)));
    else if(module_name == "div_cleaning") m_modules.push_back(std::unique_ptr<Module>(new DivCleaning(m_pd)));
    else assert(false && "Module name was not recognized");
    m_modules.back()->configureModule(in_file);
}

std::vector<std::string> ModuleHandler::getCommandLineMessages() const
{
    std::vector<std::string> messages;
    for(int i=0; i<m_modules.size(); i++){
        std::string message = m_modules[i]->commandLineMessage();
        if(!message.empty()) messages.push_back(message);
    }
    return messages;
}


void ModuleHandler::getFileOutputData(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const
{
    for(int i=0; i<m_modules.size(); i++){
        m_modules[i]->fileOutput(var_names,var_grids);
    }
}
