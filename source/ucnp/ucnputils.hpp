#ifndef UCNPUTILS_HPP
#define UCNPUTILS_HPP

#include "grid.hpp"
#include "MhdInp.hpp"
#include <filesystem>
namespace fs = std::filesystem;
#include "settings.hpp"
#include "mhd.hpp"
#include "utils.hpp"
#include "plasmadomain.hpp"
#include "antihelmholtz.hpp"

MhdInp gen_inp_grids_ucnp(const std::unique_ptr<Settings>& pms);
void genNonUniformGrids(double r_max, int Nr,std::vector<double>& dr,std::string opt);
void meshgrid(const std::vector<double>& v_x,const std::vector<double>& v_y,Grid& grid_x,Grid& grid_y);
Grid gaussian2D(Grid x,Grid y,double amp,double min,double sig_x,double sig_y,double x_cen,double y_cen);
Grid exponential2D(Grid x,Grid y,double amp,double min,double sig_x,double sig_y,double x_cen,double y_cen);

//*** PHYSICS FORMULAS - all use cgs units ***//
namespace phys
{
    double coulomb_coupling(double n,double T);
    double wigner_seitz_radius(double n);
    double plasma_freq(double n,double m);
    double einstein_freq(double n,double m);
    double nearest_coulomb_pot(double n);
    double debye_length(double n,double T);
    double screening_parameter(double n, double Te);
    double dih_temp(double n,double Te);
    double tau_exp(double sig,double m, double Te);
}

#endif