// compile: g++ -I source -I source\equationsets -I source\mhd -I source\modules -I source\modules\solar -I source\solar -I source\ucnp -I source\user-interface -std=c++20 -fopenmp -g execs\gengrids.cpp source\equationsets\*.cpp source\mhd\*.cpp source\modules\*.cpp source\modules\solar\*.cpp source\solar\*.cpp source\ucnp\*.cpp source\user-interface\*.cpp -o gengrids.exe -lm -lstdc++fs
// run: gengrids --path E:\Grant-Gorman\data-mhd\06.06.22 -s ucnp.settings -c ucnp.config --overwrite 0 --array 0
#include <iostream>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;

int main(int argc, char *argv[])
{
    std::vector<std::string> test(1);
    test[0].reserve(22);
    std::cout << "capacity before: " << test[0].capacity() << std::endl;
    test.push_back("2");
    std::cout << "capacity after: " << test[0].capacity() << std::endl;
    std::cout << "testing" << std::endl;
    return 0;
}
