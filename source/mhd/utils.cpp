#include "utils.hpp"

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

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
