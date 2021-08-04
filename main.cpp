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
#include "utils.hpp"
#include "constants.hpp"
#include "plasmadomain.hpp"

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

int main(int argc,char* argv[]){

  auto start_time = std::chrono::steady_clock::now();

  #if BENCHMARKING_ON
  Instrumentor::Get().BeginSession(out_filename);
  #endif

  std::vector<std::string> arguments = parseCommandLineArgs(argc, argv);
  std::string out_filename = arguments[0];
  std::string config_filename = arguments[1];
  std::string in_filename = arguments[2];
  int job_index = std::stoi(arguments[3]);

  PlasmaDomain simulation(out_filename.c_str(),config_filename.c_str(),job_index);

  if(in_filename.empty()){ //Default, hard coded initial condition
    simulation.hydrostaticInitialize();
    simulation.setSolarGravity(BASE_GRAV,R_SUN);
    // simulation.gaussianInitialize(1.0e-15,1.0e-12,1.0e4,3.0e4,0.05*XDIM,0.05*YDIM);
  }
  else{ //Initial condition from .state file
    simulation.readStateFile(in_filename.c_str());
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