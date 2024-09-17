#include "utils.hpp"
#include <cmath>
#include <algorithm>
#include <cassert>

//Erase all occurences of whitesapce characters in str.
//Modifies in-place.
void clearWhitespace(std::string &str)
{
  if(str.empty()) return;
  auto new_end = std::remove(str.begin(),str.end(),' ');
  for (char c : {'\t','\n','\v','\b','\r','\f','\a'}){
    new_end = std::remove(str.begin(),new_end,c);
  }
  str.erase(new_end, str.end());
}

//Wrapper for using std::getline, then removing all whitespace characters from the result
//Meant to avoid strange behavior when reading directly from a file (e.g. ifstream)
std::istream &getCleanedLine(std::istream &is, std::string &str, char delim)
{
  std::getline(is,str,delim);
  clearWhitespace(str);
  return is;
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

//Interpolates between the values of `quantity` at corresponding points (`x`,`y`)
//Evaluated at `point` (which must contain an x- and y-coordinate) using bilinear interpolation
double bilinearInterpolate(const std::vector<double> &point, const Grid &quantity, const std::vector<double> &x, const std::vector<double> &y){
  //locate grid indices in vicinity of point
  auto it_x = std::upper_bound(x.cbegin(), x.cend(), point[0]);
  auto it_y = std::upper_bound(y.cbegin(), y.cend(), point[1]);
  if(it_x == x.cend() || it_y == y.cend()) return 0.0;

  int i_1 = std::distance(x.cbegin(), it_x);
  int i_0 = std::max(i_1 - 1,0);
  int j_1 = std::distance(y.cbegin(),it_y);
  int j_0 = std::max(j_1 - 1,0);
  
  double dx0,dx1,dy0,dy1;
  dx0 = point[0] - x[i_0]; dx1 = x[i_1] - point[0];
  dy0 = point[1] - y[j_0]; dy1 = y[j_1] - point[1];

  return (quantity(i_0,j_0)*dx1*dy1 + quantity(i_1,j_0)*dx0*dy1 + quantity(i_0,j_1)*dx1*dy0 + quantity(i_1,j_1)*dx0*dy0) / ((x[i_1] - x[i_0])*(y[j_1] - y[j_0]));

}
