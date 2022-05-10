#ifndef UCNPUTILS_HPP
#define UCNPUTILS_HPP

#include "grid.hpp"
#include "MhdInp.hpp"
#include "ui_utility.hpp"
#include <filesystem>
namespace fs = std::filesystem;
#include "PlasmaSettings.hpp"
#include "settings.hpp"
#include "mhd.hpp"
#include "utils.hpp"
#include "plasmadomain.hpp"

MhdInp gen_inp_grids_ucnp(const PlasmaSettings& pms);
MhdInp gen_inp_grids_ucnp(const Settings& pms);
void genNonUniformGrids(double r_max, int Nr,std::vector<double>& dr,std::string opt);

#endif