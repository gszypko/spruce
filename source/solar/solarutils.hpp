#ifndef SOLARUTILS_HPP
#define SOLARUTILS_HPP

#include <cmath>
#include "grid.hpp"
#include "MhdInp.hpp"
#include "constants.hpp"
#include "PlasmaSettings.hpp"
#include "plasmadomain.hpp"


namespace SolarUtils {
  void SolarInitialize(Grid& rho, Grid& temp, Grid& mom_x, Grid& mom_y,
                       Grid& b_x, Grid& b_y, Grid& b_z,
                       Grid& pos_x, Grid& pos_y, Grid& grav_x, Grid& grav_y);
  
  MhdInp SolarMHDInput(const PlasmaSettings& pms);

  Grid SolarGravity(double base_gravity, double r_solar, const Grid& m_pos_y);

  //Generates potential bipolar field for component corresponding to index "index"
  //Centered s.t. origin lies at bottom middle of domain
  //Pressure scale height h, field poles at +/- l, field strength at poles b0
  Grid BipolarField(const Grid& m_pos_x, const Grid& m_pos_y, double b0, double h, int index);

  //Generates grid with hydrostatic falloff in the y-direction, with the quantity
  //"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
  //Assumes that gravity is set to vary with distance.
  Grid HydrostaticFalloff(double base_value, double scale_height, const Grid& m_pos_y);
}

#endif
