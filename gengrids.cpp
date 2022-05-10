#include <iostream>
#include "utils.hpp"
#include "settings.hpp"
#include "ucnputils.hpp"
#include "MhdInp.hpp"

namespace fs = std::filesystem;

// function notes:

int main(int argc, char *argv[])
{
    // handle input arguments
    fs::path output = getCommandLineArg(argc, argv, "-o", "--directory");
    fs::path settings = getCommandLineArg(argc, argv, "-s", "--settings");
    fs::path config = getCommandLineArg(argc, argv, "-c", "--config");
    int array = stoi(getCommandLineArg(argc, argv, "-a", "--array"));
    std::filesystem::create_directories(output);

    // load the settings file
    Settings pms(settings);
    
    // for each unique set of conditions, generate the grids and write them to file
    for (auto i : pms.get_valid_arrays()){
        fs::path set = "set_"+num2str(i+array);
        assert(!fs::exists(output/set) && "Folder path already exists - offset the array number to avoid overwriting.");
        pms.choose_array(i);
        pms.write_array_params(output/set);
        MhdInp grids = gen_inp_grids_ucnp(pms);
        grids.write_state_file(output/set);
        fs::copy_file(config,output/set/config);
    }
    
    std::cout << "complete" << std::endl;
    return 0;
}
