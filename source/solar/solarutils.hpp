#ifndef SOLARUTILS_HPP
#define SOLARUTILS_HPP

#include "grid.hpp"

namespace SolarUtils {

  Grid SolarGravity(double base_gravity, double r_solar, const Grid& m_pos_y);

  //Generates potential bipolar field for component corresponding to index "index"
  //Centered s.t. origin lies at bottom middle of domain
  //Pressure scale height h, field poles at +/- l, field strength at poles b0
  Grid BipolarField(const Grid& m_pos_x, const Grid& m_pos_y, double b0, double h, int index);

  //Generates grid with hydrostatic falloff in the y-direction, with the quantity
  //"base_value" at y=0. Assumes isothermal atmosphere with temperature "iso_temp".
  //Assumes that gravity is set to vary with distance.
  Grid HydrostaticFalloff(double base_value, double scale_height, const Grid& m_pos_y);

  //Generates gaussian initial condition for a variable, centered at middle of grid
  //std_dev_x and std_dev_y are the standard deviation of the distribution in the x
  //and y directions, in units of grid cell widths
  Grid GaussianGrid(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y, double center_x, double center_y);
  Grid GaussianGridRotated(int xdim, int ydim, double min, double max, double std_dev_x, double std_dev_y, double center_x, double center_y, double angle_deg);
  //In version where distribution center is unspecified, placed at center of domain
  Grid GaussianGrid(const Grid& pos_x, const Grid& pos_y, double min, double max, double std_dev_x, double std_dev_y);
}

#endif
