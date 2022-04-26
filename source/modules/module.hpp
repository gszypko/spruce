//module.hpp
//Defines abstract base class Module
//All specific Modules must be defined/implemented as
//a derived class of Module

#ifndef MODULE_HPP
#define MODULE_HPP

#include "constants.hpp"
#include "utils.hpp"
#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

class PlasmaDomain;

class Module {
    public:
        Module(PlasmaDomain &pd);
        virtual ~Module() {}
        void configureModule(std::ifstream &in_file);
        //Evolve the system in time by dt according to the Module
        //Note that m_pd.propagateChanges() is NOT automatically applied;
        //this must be done manually such that all variables remains consistent
        virtual void iterateModule(double dt);
        //Perform any processing that occurs before the main iteration step
        //All modules are guaranteed to run this function before any of them
        //run iterateModule() or postIterateModule()
        virtual void preIterateModule(double dt);
        //Perform any processing that occurs after the main iteration step
        //and after the evolution of the main MHD equations by m_pd
        //All modules are guaranteed to run this function after all of them
        //run iterateModule() and postIterateModule()
        virtual void postIterateModule(double dt);
        //Returns a short message to print to stdout related to the most recent iteration
        //of the Module. Should not include any line breaks.
        //Default behavior is no message; override in derived Module classes to customize.
        virtual std::string commandLineMessage() const;
        //Returns data to write to the PlasmaDomain output file corresponding to
        //the most recent iteration of the Module.
        //Any variables to output must have their names appended to var_names
        //and the corresponding Grids appended to var_grids, in the same order. 
        //Default behavior is no data; override in derived Module classes to customize.
        virtual void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const;
    protected:
        PlasmaDomain& m_pd;
        //Apply the values (in rhs) to the appropriate Module configs (in lhs)
        virtual void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) = 0;
};

#endif