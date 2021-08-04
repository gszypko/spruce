#include "constants.hpp"
#include "utils.hpp"
#include "grid.hpp"
#include <omp.h>
#include <cmath>
#include <vector>
#include <limits>
#include <string>
#include <algorithm>
#include <sstream>

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

//Generates gaussian initial condition for a variable, centered at middle of grid
//std_dev_x and std_dev_y are the standard deviation of the distribution in the x
//and y directions, in units of grid cell widths
Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  std::vector<double> gauss_x(xdim), gauss_y(ydim);
  // double sigmax = 0.05*xdim;
  // double sigmay = 0.05*ydim;
  for(int i=0; i<xdim; i++){
    gauss_x[i] = std::exp(-0.5*std::pow(((double)i-0.5*(double)(xdim-1))/std_dev_x,2.0));
  }
  for(int j=0; j<ydim; j++){
    gauss_y[j] = std::exp(-0.5*std::pow(((double)j-0.5*(double)(ydim-1))/std_dev_y,2.0));
  }
  Grid result(xdim,ydim);
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      result(i,j) = gauss_x[i]*gauss_y[j];
    }
  }
  return result;
}

//Generates potential bipolar field for component corresponding to index "index"
//Centered s.t. origin lies at bottom middle of domain
//Pressure scale height h, field poles at +/- l, field strength at poles b0
Grid BipolarField(int xdim, int ydim, double b0, double h, double dx, double dy, int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result = Grid::Zero(xdim, ydim);
  #pragma omp parallel for collapse(2)
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      double x = (i - (double)(xdim-1)*0.5)*dx;
      double y = j*dy;
      if(index == 0) result(i,j) = b0*std::exp(-0.5*y/h)*std::cos(0.5*x/h);
      else result(i,j) = -b0*std::exp(-0.5*y/h)*std::sin(0.5*x/h);
    }
  }
  return result;
}

//Generates grid with hydrostatic falloff in the y-direction, with the quantity
//"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
//Assumes that gravity is set to vary with distance.
Grid HydrostaticFalloff(double base_value, double scale_height, int xdim, int ydim, double dy){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result = Grid::Zero(xdim, ydim);
  #pragma omp parallel for collapse(2)
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      double y = j*dy;
      //result(i,j) = base_value*std::exp(-y/scale_height);
      result(i,j) = base_value*std::exp(M_SUN*GRAV_CONST/BASE_GRAV/scale_height*(1/(y+R_SUN) - 1/R_SUN));
    }
  }
  return result;
}

//Erase all occurences of ' ', '\t', and '\n' in str.
//Modifies in-place.
void clearWhitespace(std::string &str)
{
  if(str.empty()) return;
  auto new_end = std::remove(str.begin(),str.end(),' ');
  new_end = std::remove(str.begin(),new_end,'\t');
  new_end = std::remove(str.begin(),new_end,'\n');
  str.erase(new_end, str.end());
}

//Splits string into vector of substrings, delimited by delim
std::vector<std::string> splitString(const std::string &str, const char delim)
{
  std::vector<std::string> result;
  std::istringstream ss(str);
  std::string element;
  while(std::getline(ss,element,delim)){
    if(!element.empty()) result.push_back(element);
  }
  return result;
}

std::vector<std::string> parseCommandLineArgs(int argc, char* argv[])
{
  std::vector<std::string> arguments(argv+1, argv+argc);
  std::vector<std::string> result(3);
  int num_args = arguments.size();
  assert(num_args <= 2*result.size() && "No more than one argument per flag");
  for(int i=0; i<num_args; i++){
    assert(arguments[i][0] == '-' && i+1<num_args && "Arguments must be given as flag followed by non-flag");
    if(arguments[i] == "-o" || arguments[i] == "--output") result[0] = arguments[i+1];
    else if(arguments[i] == "-c" || arguments[i] == "--config") result[1] = arguments[i+1];
    else if(arguments[i] == "-s" || arguments[i] == "--state") result[2] = arguments[i+1];
    i++;
  }
  if(result[0].empty()) result[0] = "output";
  if(result[1].empty()) result[1] = "default.config";
  return result;
}
