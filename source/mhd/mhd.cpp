#include "grid.hpp"
#include "plasmadomain.hpp"

void mhd(const Grid& rho, const Grid& temp, const Grid& mom_x, const Grid& mom_y,
         const Grid& b_x, const Grid& b_y, const Grid& b_z, const Grid& pos_x, const Grid& pos_y,
         const Grid& grav_x, const Grid& grav_y, const char* output_pathname, const char* config_filename)
{
  PlasmaDomain simulation(output_pathname,config_filename);
  simulation.initialize(rho, temp, mom_x, mom_y, b_x, b_y, b_z, pos_x, pos_y, grav_x, grav_y);
 simulation.run();
}

void mhd(const char* state_filename, const char* output_pathname, const char* config_filename)
{
  PlasmaDomain simulation(output_pathname,config_filename);
  simulation.readStateFile(state_filename);
  simulation.run();
}