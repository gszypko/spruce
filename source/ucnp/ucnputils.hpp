#ifndef UCNPUTILS_HPP
#define UCNPUTILS_HPP


#include "MhdInp.hpp"
#include "ui_utility.hpp"
#include <filesystem>
#include "PlasmaSettings.hpp"
#include "mhd.hpp"
#include "utils.hpp"
#include "plasmadomain.hpp"

MhdInp gen_inp_grids_ucnp(const PlasmaSettings& pms);

#endif