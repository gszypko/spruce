#ifndef MHD_HPP
#define MHD_HPP

#include "grid.hpp"

//Explicit variable specification
void mhd(const Grid& rho, const Grid& temp, const Grid& mom_x, const Grid& mom_y,
         const Grid& b_x, const Grid& b_y, const Grid& b_z, const Grid& pos_x, const Grid& pos_y,
         const Grid& grav_x, const Grid& grav_y, const char* output_pathname, const char* config_filename);

//Variable specification from .state file
void mhd(const char* state_filename, const char* output_pathname, const char* config_filename);

#endif