#include "utils.hpp"

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

std::string getCommandLineArg(int argc, char* argv[], std::string short_flag, std::string long_flag)
{
  std::vector<std::string> arguments(argv+1, argv+argc);
  std::string result;
  int num_args = argc - 1;
  for(int i=0; i<num_args; i++){
    assert(arguments[i][0] == '-' && i+1<num_args && "Arguments must be given as flag followed by non-flag");
    if(arguments[i] == short_flag || arguments[i] == long_flag){
      assert(result.empty() && "Command line flags cannot be used more than once");
      result = arguments[i+1];
    }
    i++;
  }
  return result;
}
