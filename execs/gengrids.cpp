// compile: g++ -I source -I source\mhd -I source\modules\solar -I source\modules -I source\solar -I source\ucnp -I source\user-interface -std=c++20 -fopenmp -O3 execs\gengrids.cpp source\mhd\*.cpp source\modules\*.cpp source\modules\solar\*.cpp source\solar\*.cpp source\ucnp\*.cpp source\user-interface\*.cpp -o execs\gengrids.exe -lm -lstdc++fs
// run: execs\gengrids -o E:\Grant-Gorman\data-mhd\05.12.22 -s ucnp.settings -c ucnp.config -f 0 -a 0
#include <iostream>
#include "utils.hpp"
#include "settings.hpp"
#include "ucnputils.hpp"
#include "MhdInp.hpp"

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    // handle input arguments
    fs::path output = getCommandLineArg(argc, argv, "-o", "--directory");
    fs::path settings = getCommandLineArg(argc, argv, "-s", "--settings");
    fs::path config = getCommandLineArg(argc, argv, "-c", "--config");
    int array = stoi(getCommandLineArg(argc, argv, "-a", "--array"));
    std::string overwrite_str = getCommandLineArg(argc, argv, "-f", "--flag");

    assert(settings.extension().string()==".settings" && fs::exists(settings) && "Error: settings file must exist and have extension .settings");
    assert(config.extension().string()==".config" && fs::exists(config) && "Error: config file must exist and have extension .config");

    int overwrite{0};
    if (!overwrite_str.empty()) overwrite = stoi(overwrite_str);
    assert((overwrite==0 || overwrite==1) && "Error: the overwrite flag is logical and must be 0 or 1");

    std::filesystem::create_directories(output);

    // load the settings file
    Settings pms(settings);
    
    // for each unique set of conditions, generate the grids and write them to file
    for (auto i : pms.get_valid_arrays()){
        fs::path set_path = output/("set_"+num2str(i+array));
        if (overwrite==0) assert(!fs::exists(set_path) && "Error: folder already exists and overwrite_flag=0");
        
        pms.choose_array(i);
        pms.write_array_params(set_path);
        
        MhdInp grids = gen_inp_grids_ucnp(pms);
        grids.write_state_file(set_path);

        fs::path new_config_path = set_path/config;
        if (new_config_path.extension().string()==".config" && fs::exists(new_config_path)) fs::remove(new_config_path);
        fs::copy_file(config,new_config_path);
    }
    
    std::cout << "Complete." << std::endl;
    return 0;
}
