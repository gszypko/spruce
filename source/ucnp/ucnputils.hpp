#ifndef UCNPUTILS_HPP
#define UCNPUTILS_HPP

#include "grid.hpp"
#include "MhdInp.hpp"
#include "ui_utility.hpp"
#include <filesystem>
#include "PlasmaSettings.hpp"
#include "mhd.hpp"
#include "utils.hpp"
#include "plasmadomain.hpp"

MhdInp gen_inp_grids_ucnp(const PlasmaSettings& pms);

void genNonUniformGrids(double r_max, double dr_min, int Nr,std::vector<double>& dr);

#endif