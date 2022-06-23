// compile: g++ -I source -I source\equationsets -I source\mhd -I source\modules -I source\modules\solar -I source\solar -I source\ucnp -I source\user-interface -std=c++20 -fopenmp -g execs\gengrids.cpp source\equationsets\*.cpp source\mhd\*.cpp source\modules\*.cpp source\modules\solar\*.cpp source\solar\*.cpp source\ucnp\*.cpp source\user-interface\*.cpp -o gengrids.exe -lm -lstdc++fs
// run: gengrids --path E:\Grant-Gorman\data-mhd\06.06.22 -s ucnp.settings -c ucnp.config --overwrite 0 --array 0
#include <iostream>
#include "utils.hpp"
#include "settings.hpp"
#include "ucnp.hpp"
#include "ucnputils.hpp"
#include "MhdInp.hpp"

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    std::unique_ptr<Settings> ptr = Settings::spawn_settings("ucnp","ucnp.settings");
    return 0;
}
