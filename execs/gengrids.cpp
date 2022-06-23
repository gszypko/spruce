// compile: g++ -I source -I source\equationsets -I source\mhd -I source\modules -I source\modules\solar -I source\solar -I source\ucnp -I source\user-interface -std=c++20 -fopenmp -g execs\gengrids.cpp source\equationsets\*.cpp source\mhd\*.cpp source\modules\*.cpp source\modules\solar\*.cpp source\solar\*.cpp source\ucnp\*.cpp source\user-interface\*.cpp -o gengrids.exe -lm -lstdc++fs
// run: gengrids --path E:\Grant-Gorman\data-mhd\06.06.22 -s ucnp.settings -c ucnp.config --overwrite 0 --array 0
#include <iostream>
#include "utils.hpp"
#include "settings.hpp"
#include "ucnputils.hpp"
#include "MhdInp.hpp"

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    // read in command line arguments
    fs::path output = getCommandLineArg(argc, argv, "-p", "--path");
    fs::path settings = getCommandLineArg(argc, argv, "-s", "--settings");
    fs::path config = getCommandLineArg(argc, argv, "-c", "--config");
    std::string array_str = getCommandLineArg(argc, argv, "-a", "--array");
    std::string overwrite_str = getCommandLineArg(argc, argv, "-o", "--overwrite");
    // ensure that .settings and .config files exist and have the correct extensions
    assert(settings.extension().string()==".settings" && fs::exists(settings) && "Error: settings file must exist and have extension .settings");
    assert(config.extension().string()==".config" && fs::exists(config) && "Error: config file must exist and have extension .config");
    // handle command line inputs
    int overwrite{0},array{0};
    if (!overwrite_str.empty()) overwrite = stoi(overwrite_str);
    if (!array_str.empty()) array = stoi(array_str);
    assert((overwrite==0 || overwrite==1) && "Error: the overwrite flag is logical and must be 0 or 1");
    // create output directory structure
    std::filesystem::create_directories(output);
    // load the settings file
    std::unique_ptr<Settings> pms = Settings::spawn_settings("ucnp",settings);
    pms->print_task_array_range();
    // for each unique set of conditions, generate the grids and write them to file
    for (auto i : pms->task_array()){
        // choose array and write array parameters to output directory
        pms->choose_array(i);
        fs::path set_path = output/pms->set_path(array);
        if (overwrite==0) assert(!fs::exists(set_path) && "Error: folder already exists and overwrite_flag=0");
        pms->write_array_params(set_path,"plasma");
        // create grids and write to .state file
        MhdInp grids = gen_inp_grids_ucnp(pms);
        grids.write_state_file(set_path);
        // handle .config file
        grids.read_config(config);
        for (const auto& name : pms->names()){
            if (grids.is_config(name)){
                grids.update_config(name,Settings::num2str(pms->getvar(name)));
            }
        }
        grids.write_config_file(set_path);
    }    
    return 0;
}
