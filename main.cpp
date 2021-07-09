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

void sigfpe_handler( int signal_num ){
  std::cout << "Floating point exception (" << signal_num << "), halting execution.\n";
  exit(signal_num);
}

int main(int argc,char* argv[]){
  std::signal(SIGFPE, sigfpe_handler);

  std::string out_filename;
  if(argc == 1) out_filename = "output";
  else out_filename = std::string(argv[1]);

  #if BENCHMARKING_ON
  Instrumentor::Get().BeginSession(out_filename);
  #endif

  auto start_time = std::chrono::steady_clock::now();

  //Set up initial conditions
  PlasmaDomain simulation(XDIM,YDIM,out_filename.c_str());
  double time_in;
  if(argc == 3){
    const char* in_filename = argv[2];
    simulation.readStateFile(in_filename);
  }
  else{ //Hard coded initial condition
    simulation.hydrostaticInitialize();
    double b0 = B0; //magnetic field strength at base in G
    // double iso_temp = TEMP_CHROMOSPHERE; //isothermal initial temp in K
    // b_x = BipolarField(XDIM, YDIM, b0, scale_height, 0);
    // b_y = BipolarField(XDIM, YDIM, b0, scale_height, 1);
    // b_z = Grid::Zero(XDIM,YDIM);

    // //Simple Gaussian test case
    // rho = GaussianGrid(XDIM, YDIM, 1.0, 5.0);
    // mom_x = Grid::Zero(XDIM,YDIM); //x momentum density
    // mom_y = Grid::Zero(XDIM,YDIM); //y momentum density
    // temp = GaussianGrid(XDIM, YDIM, 1.0, 2.0); //temperature
    // double b0 = 0.1;
    // double h = (XDIM*DX)/(4.0*PI);
    // b_x = BipolarField(XDIM, YDIM, b0, h, 0);
    // b_y = BipolarField(XDIM, YDIM, b0, h, 1);
    // b_z = Grid::Zero(XDIM,YDIM);

    //Isothermal hydrostatic initial condition
    // rho = HydrostaticFalloff(base_rho,scale_height,XDIM,YDIM);
    // mom_x = Grid::Zero(XDIM,YDIM); //x momentum density
    // mom_y = Grid::Zero(XDIM,YDIM); //y momentum density
    // temp = iso_temp*Grid::Ones(XDIM,YDIM); //temperature

    time_in = 0.0; //default initial time is zero
  }

  simulation.setSolarGravity(BASE_GRAV,R_SUN);
  simulation.recomputeDerivedVariables();

  //Outputs domain dimensions and magnetic field
  simulation.outputPreamble();
  simulation.outputCurrentState();

  for (int iter = 0; iter < NT; iter++){
    #if BENCHMARKING_ON
    InstrumentationTimer timer((std::string("iteration ") + std::to_string(iter)).c_str());
    #endif

    // Enforce rigid boundaries
    // if(YBOUND1 == WALL || YBOUND2 == WALL){
    //   if(YBOUND1 == WALL) for(int i=0; i<XDIM; i++){
    //     mom_x(i,0) = 0.0;
    //     mom_y(i,0) = 0.0;
    //   }
    //   if(YBOUND2 == WALL) for(int i=0; i<XDIM; i++){
    //     mom_x(i,YDIM-1) = 0.0;
    //     mom_y(i,YDIM-1) = 0.0;
    //   }
    // }
    // if(XBOUND1 == WALL || XBOUND2 == WALL){
    //   if(XBOUND1 == WALL) for(int j=0; j<YDIM; j++){
    //     mom_x(0,j) = 0.0;
    //     mom_y(0,j) = 0.0;
    //   }
    //   if(XBOUND2 == WALL) for(int j=0; j<YDIM; j++){
    //     mom_x(XDIM-1,j) = 0.0;
    //     mom_y(XDIM-1,j) = 0.0;
    //   }
    // }

    simulation.advanceTime();

    //Clamping wall boundary values
    // if(YBOUND1 == WALL || YBOUND2 == WALL){
    //   if(YBOUND1 == WALL) for(int i=0; i<XDIM; i++){
    //     mom_x_next(i,0) = 0.0;
    //     mom_y_next(i,0) = 0.0;
    //     rho_next(i,0) = rho(i,0);
    //     energy_next(i,0) = energy(i,0);
    //   }
    //   if(YBOUND2 == WALL) for(int i=0; i<XDIM; i++){
    //     mom_x_next(i,YDIM-1) = 0.0;
    //     mom_y_next(i,YDIM-1) = 0.0;
    //     rho_next(i,YDIM-1) = rho(i,YDIM-1);
    //     energy_next(i,YDIM-1) = energy(i,YDIM-1);
    //   }
    // }
    // if(XBOUND1 == WALL || XBOUND2 == WALL){
    //   if(XBOUND1 == WALL) for(int j=0; j<YDIM; j++){
    //     mom_x_next(0,j) = 0.0;
    //     mom_y_next(0,j) = 0.0;
    //     rho_next(0,j) = rho(0,j);
    //     energy_next(0,j) = energy(0,j);
    //   }
    //   if(XBOUND2 == WALL) for(int j=0; j<YDIM; j++){
    //     mom_x_next(XDIM-1,j) = 0.0;
    //     mom_y_next(XDIM-1,j) = 0.0;
    //     rho_next(XDIM-1,j) = rho(XDIM-1,j);
    //     energy_next(XDIM-1,j) = energy(XDIM-1,j);
    //   }
    // }

    simulation.catchUnderdensity();
    simulation.recomputeTemperature();
    simulation.recomputeDerivedVariables();

    if(iter%OUTPUT_INTERVAL == 0){
      simulation.outputCurrentState();
    }

    simulation.writeStateFile();
  }

  //Clean up state files
  rename((out_filename+std::to_string((NT-1)%2)+".state").c_str(),
          (out_filename+".state").c_str());
  remove((out_filename+std::to_string((NT-2)%2)+".state").c_str());

  //Information about run time
  std::cout << "\rIterations: " << NT << "/" << NT << "\n";
  auto stop_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop_time - start_time);
  double minutes = (int)(duration.count()/60.0);
  double seconds = duration.count() - 60*minutes;
  // out_file << "runtime=" << minutes << "min" << seconds << "sec";
  std::cout << "\rTotal runtime: " << minutes << " min " << seconds << " sec (approx. " 
    << (double)duration.count()/(double)NT << " sec per iteration)" << std::endl;
  // out_file.close();

  #if BENCHMARKING_ON
  Instrumentor::Get().EndSession();
  #endif
}