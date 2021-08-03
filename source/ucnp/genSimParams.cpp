#include "ucnp_ui.hpp"

int main(int argc, char *argv[])
{
    std::string settings_path = "source/ucnp/ucnp.settings";
    std::vector<std::vector<std::string>> test = readCSV(settings_path);

    enum Vars {v_n,v_sigX,v_sigY,v_Ti,v_Te,v_dBdx,v_tmax,v_xMax,v_yMax,v_Nx,v_Ny,v_mi,v_num};
    std::vector<std::vector<double>> vars(test.size());

    for (int i = 0; i < test.size(); i++){
        for (int j = 0; j < test[i].size(); j++){
            vars[i].push_back(atof(test[i][j].c_str()));
        }
    }


    return 0;
}