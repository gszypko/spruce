#ifndef DERIVS_HPP
#define DERIVS_HPP

#include "grid.hpp"

//Compute surface values from cell-centered values using Barton's method
//Meant to be used for transport terms only
//Result indexed s.t. element i,j indicates surface between i,j and i-1,j
//if "index"==0, or i,j and i,j-1 if "index"==1
Grid upwind_surface(const Grid &cell_center, const Grid &vel, const int index);

Grid transport_derivative1D(const Grid &quantity, const Grid &vel, const int index);

//Compute single-direction divergence term for non-transport term (central differencing)
Grid derivative1D(const Grid &quantity, const int index);

//Compute divergence term for simulation parameter "quantity"
//"quantity","vx","vy" used for transport term
//Non-transport terms contained in "nontransp_x", "nontransp_y"
Grid divergence(const Grid &quantity, const Grid &nontransp_x, const Grid &nontransp_y, const Grid &vx, const Grid &vy);

//Compute single-direction second derivative
Grid second_derivative1D(const Grid &quantity, const int index);

//Computes Laplacian (del squared) of "quantity"
Grid laplacian(const Grid &quantity);

#endif
