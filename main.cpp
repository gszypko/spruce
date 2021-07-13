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

  //Set up initial conditions
  PlasmaDomain simulation(XDIM,YDIM,DX,DY,out_filename.c_str()); //Note: dimensions given here can be overridden by state file
  simulation.setDefaultSettings(); //Applies #define macros to member variables of PlasmaDomain object

  if(argc == 3){ //Initial condition from .state file
    const char* in_filename = argv[2];
    simulation.readStateFile(in_filename);
  }
  else{ //Default, hard coded initial condition
    simulation.hydrostaticInitialize();
    // simulation.gaussianInitialize();
  }

  simulation.setSolarGravity(BASE_GRAV,R_SUN);

  simulation.outputPreamble();
  simulation.outputCurrentState();

  for (int iter = 0; iter < NT; iter++){
    #if BENCHMARKING_ON
    InstrumentationTimer timer((std::string("iteration ") + std::to_string(iter)).c_str());
    #endif

    simulation.advanceTime();

    if(iter%OUTPUT_INTERVAL == 0){
      simulation.outputCurrentState();
    }

    simulation.writeStateFile();
  }

  simulation.cleanUpStateFiles();

  //Information about run time
  std::cout << "\rIterations: " << NT << "/" << NT << "\n";
  auto stop_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time);
  double minutes = (int)(duration.count()/60.0);
  double seconds = duration.count() - 60*minutes;
  std::cout << "\rTotal runtime: " << minutes << " min " << seconds << " sec (approx. " 
    << (double)duration.count()/(double)NT << " sec per iteration)" << std::endl;

  #if BENCHMARKING_ON
  Instrumentor::Get().EndSession();
  #endif
}