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
#include "mhd.hpp"
#include "solarutils.hpp"

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

  if(in_filename.empty()){ //Default, hard coded initial condition
    Grid rho, temp, mom_x, mom_y, b_x, b_y, b_z, pos_x, pos_y, grav_x, grav_y;
    SolarUtils::SolarInitialize(rho, temp, mom_x, mom_y, b_x, b_y, b_z, pos_x, pos_y, grav_x, grav_y);
    mhd(rho, temp, mom_x, mom_y, b_x, b_y, b_z, pos_x, pos_y, grav_x, grav_y, out_filename.c_str(), config_filename.c_str());
  }
  else{ //Initial condition from .state file
    mhd(in_filename.c_str(),out_filename.c_str(),config_filename.c_str());
  }

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