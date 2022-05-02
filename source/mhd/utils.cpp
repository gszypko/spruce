#include "utils.hpp"

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

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
std::istream &getCleanedLine(std::istream &is, std::string &str)
{
  std::getline(is,str);
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
