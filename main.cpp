//mhdtoy
//Greg Szypko

#include <iostream>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <csignal>
#include <cstdio>
#include <cstring>
#include <chrono>
#include <string>
#include <vector>
#include "grid.hpp"
#include "derivs.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "plasmadomain.hpp"

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

int main(int argc,char* argv[]){

  std::string out_filename;
  if(argc == 1) out_filename = "output";
  else out_filename = std::string(argv[1]);

  #if BENCHMARKING_ON
  Instrumentor::Get().BeginSession(out_filename);
  #endif

  auto start_time = std::chrono::steady_clock::now();

  //Initialize PlasmaDomain with settings from define macros
  PlasmaDomain simulation(XDIM,YDIM,DX,DY,out_filename.c_str()); //Note: dimensions given here can be overridden by state file
  simulation.setDefaultSettings(); //Applies #define macros to member variables of PlasmaDomain object

  // //Initialize PlasmaDomain with settings from .settings file (under construction)
  // PlasmaDomain simulation(out_filename.c_str(),"default.settings");

  if(argc == 3){ //Initial condition from .state file
    const char* in_filename = argv[2];
    simulation.readStateFile(in_filename);
  }
  else{ //Default, hard coded initial condition
    simulation.hydrostaticInitialize();
    simulation.setSolarGravity(BASE_GRAV,R_SUN);
    // simulation.gaussianInitialize(1.0e-15,1.0e-12,1.0e4,3.0e4,0.05*XDIM,0.05*YDIM);
  }

  simulation.run();

  //Information about run time
  auto stop_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time);
  double minutes = (int)(duration.count()/60.0);
  double seconds = duration.count() - 60*minutes;
  std::cout << "\rTotal runtime: " << minutes << " min " << seconds << " sec" << std::endl;

  #if BENCHMARKING_ON
  Instrumentor::Get().EndSession();
  #endif
}