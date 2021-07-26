//derivs.cpp
//Defines differentiation functions for PlasmaDomain

#include "constants.hpp"
#include "plasmadomain.hpp"
#include "grid.hpp"
#include <string>
#include <cmath>
#include <omp.h>

#if BENCHMARKING_ON
#include "instrumentor.hpp"
#endif

//Compute surface values from cell-centered values using Barton's method
//Meant to be used for transport terms only
//Result indexed s.t. element i,j indicates surface between i,j and i-1,j
//if "index"==0, or i,j and i,j-1 if "index"==1
Grid PlasmaDomain::upwindSurface(const Grid &cell_center, const Grid &vel, const int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  int xdim = m_xdim+1-index;
  int ydim = m_ydim+index;
  Grid cell_surface = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #if BENCHMARKING_ON
    InstrumentationTimer timer("loop thread");
    #endif
    #pragma omp for collapse(2)
    for (int i = 0; i < xdim; i++){
      for(int j = 0; j < ydim; j++){
        //Handle direction of cell_center being considered (i.e. index for differencing)
        int i2 = i, j2 = j;
        int i0, i1, i3, j0, j1, j3;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j2; j1 = j2; j3 = j2;
          i0 = i2-2; i1 = i2-1; i3 = i2+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            if(i2 == m_xdim) continue;
            //Here, explicitly need macro m_xdim instead of xdim
            i0 = (i0+m_xdim)%m_xdim;
            i1 = (i1+m_xdim)%m_xdim;
            i3 = (i3+m_xdim)%m_xdim;
          }
          if(x_bound_1 == BoundaryCondition::Wall){
            if(i2 == 0){
              cell_surface(i2,j2) = 0.0;
              continue;
            } else if(i2 == 1){
              i0 = i1;
            }
          } else if(x_bound_1 == BoundaryCondition::Open){
            if(i2 == 0){
              cell_surface(i2,j2) = 1.5*cell_center(i2,j2) - 0.5*cell_center(i3,j3); //lerp
              continue;
            } else if(i2 == 1){
              i0 = i1;
            }
          }
          if(x_bound_2 == BoundaryCondition::Wall){
            if(i2 == m_xdim){
              cell_surface(i2,j2) = 0.0;
              continue;
            } else if(i2 == m_xdim - 1){
              i3 = i2;
            }
          } else if(x_bound_2 == BoundaryCondition::Open){
            if(i2 == m_xdim){
              cell_surface(i2,j2) = 1.5*cell_center(i1,j1) - 0.5*cell_center(i0,j0); //lerp
              continue;
            } else if(i2 == m_xdim - 1){
              i3 = i2;
            }
          }
        }
        else{
          //Handle Y boundary conditions
          i0 = i2; i1 = i2; i3 = i2;
          j0 = j2-2; j1 = j2-1; j3 = j2+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            if(j2 == m_ydim) continue;
            //Here, explicitly need macro m_ydim instead of ydim
            j0 = (j0+m_ydim)%m_ydim;
            j1 = (j1+m_ydim)%m_ydim;
            j3 = (j3+m_ydim)%m_ydim;
          }
          if(y_bound_1 == BoundaryCondition::Wall){
            if(j2 == 0){
              cell_surface(i2,j2) = 0.0;
              continue;
            } else if(j2 == 1){
              j0 = j1;
            }
          } else if(y_bound_1 == BoundaryCondition::Open){
            if(j2 == 0){
              cell_surface(i2,j2) = 1.5*cell_center(i2,j2) - 0.5*cell_center(i3,j3); //lerp
              continue;
            } else if(j2 == 1){
              j0 = j1;
            }
          }
          if(y_bound_2 == BoundaryCondition::Wall){
            if(j2 == m_ydim){
              cell_surface(i2,j2) = 0.0;
              continue;
            } else if(j2 == m_ydim - 1){
              j3 = j2;
            }
          } else if(y_bound_2 == BoundaryCondition::Open){
            if(j2 == m_ydim){
              cell_surface(i2,j2) = 1.5*cell_center(i1,j1) - 0.5*cell_center(i0,j0); //lerp
              continue;
            } else if(j2 == m_ydim - 1){
              j3 = j2;
            }
          }
        }
        //Apply Barton's method
        double d1, d2, d3;
        d2 = 0.5*(cell_center(i1,j1)+cell_center(i2,j2));
        if(0.5*(vel(i1,j1)+vel(i2,j2))>=0.0){
          d3 = cell_center(i1,j1);
          d1 = 1.5*cell_center(i1,j1) - 0.5*cell_center(i0,j0);
          if(cell_center(i2,j2) <= cell_center(i1,j1)){
            cell_surface(i2,j2) = std::min(d3,std::max(d1,d2));
          } else { //cell_center(i2,j2) > cell_center(i1,j1)
            cell_surface(i2,j2) = std::max(d3,std::min(d1,d2));
          }
        } else { //vel(i1,j1)<0.0
          d3 = cell_center(i2,j2);
          d1 = 1.5*cell_center(i2,j2) - 0.5*cell_center(i3,j3);
          if(cell_center(i2,j2) <= cell_center(i1,j1)){
            cell_surface(i2,j2) = std::max(d3,std::min(d1,d2));
          } else { //cell_center(i2,j2) > cell_center(i1,j1)
            cell_surface(i2,j2) = std::min(d3,std::max(d1,d2));
          }
        }
      }
    }
  }
  return cell_surface;
}

Grid PlasmaDomain::transportDerivative1D(const Grid &quantity, const Grid &vel, const int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid surf_quantity = upwindSurface(quantity, vel, index);
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  double denom = m_dx*(1-index) + m_dy*(index);
  Grid div = Grid::Ones(xdim,ydim);
  #pragma omp parallel
  {
    #if BENCHMARKING_ON
    InstrumentationTimer timer("loop thread");
    #endif
    #pragma omp for collapse(2)
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        int i0, i1, i2, i2surf, j0, j1, j2, j2surf; //Need separate indices for surface and vel
        i1 = i; j1 = j;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j1; j2 = j1; j2surf = j1;
          i0 = i1-1; i2 = i1+1; i2surf = i1+1;
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+m_xdim)%m_xdim;
            i2 = (i2+m_xdim)%m_xdim;
            i2surf = (i2surf+m_xdim)%m_xdim;
          }
          if(x_bound_1 == BoundaryCondition::Wall){
            // if(i1 == 0 || i1 == 1){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(i1 == 0){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(x_bound_1 == BoundaryCondition::Open){
            if(i1 == 0) i0 = i1;
          }
          if(x_bound_2 == BoundaryCondition::Wall){
            // if(i1 == m_xdim-1 || i1 == m_xdim-2){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(i1 == m_xdim-1){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(x_bound_2 == BoundaryCondition::Open){
            if(i1 == m_xdim-1) i2 = i1;
          }
          // if(x_bound_1 == BoundaryCondition::Wall || x_bound_1 == BoundaryCondition::Open){
          //   if(i1 == 0) i0 = i1;
          // }
          // if(x_bound_2 == BoundaryCondition::Wall || x_bound_2 == BoundaryCondition::Open){
          //   if(i1 == m_xdim-1) i2 = i1;
          // }
        }
        else{
          //Handle Y boundary conditions
          i0 = i1; i2 = i1; i2surf = i1;
          j0 = j1-1; j2 = j1+1; j2surf = j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+m_ydim)%m_ydim;
            j2 = (j2+m_ydim)%m_ydim;
            j2surf = (j2surf+m_ydim)%m_ydim;
          }
          if(y_bound_1 == BoundaryCondition::Wall){
            // if(j1 == 0 || j1 == 1){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(j1 == 0){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(y_bound_1 == BoundaryCondition::Open){
            if(j1 == 0) j0 = j1;
          }
          if(y_bound_2 == BoundaryCondition::Wall){
            // if(j1 == m_ydim-1 || j1 == m_ydim-2){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(j1 == m_ydim-1){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(y_bound_2 == BoundaryCondition::Open){
            if(j1 == m_ydim-1) j2 = j1;
          }
          // if(y_bound_1 == BoundaryCondition::Wall || y_bound_1 == BoundaryCondition::Open){
          //   if(j1 == 0) j0 = j1;
          // }
          // if(y_bound_2 == BoundaryCondition::Wall || y_bound_2 == BoundaryCondition::Open){
          //   if(j1 == m_ydim-1) j2 = j1;
          // }
        }
        div(i1,j1) = (surf_quantity(i2surf,j2surf)*0.5*(vel(i1,j1)+vel(i2,j2))
                  - surf_quantity(i1,j1)*0.5*(vel(i1,j1)+vel(i0,j0)))/denom;
      }
    }
  }
  return div;
}

//Compute single-direction divergence term for non-transport term (central differencing)
Grid PlasmaDomain::derivative1D(const Grid &quantity, const int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Zero(xdim,ydim);
  double denom = 2.0*(m_dx*(1-index) + m_dy*(index));
  #pragma omp parallel
  {
    #if BENCHMARKING_ON
    InstrumentationTimer timer("loop thread");
    #endif
    #pragma omp for collapse(2)
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        int i0, i1, i2, j0, j1, j2;
        i1 = i; j1 = j;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j1; j2 = j1;
          i0 = i1-1; i2 = i1+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+xdim)%xdim;
            i2 = (i2+xdim)%xdim;
          }
          if(x_bound_1 == BoundaryCondition::Wall){
            // if(i1 == 0 || i1 == 1){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(i1 == 0){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(x_bound_1 == BoundaryCondition::Open){
            if(i1 == 0) i0 = i1;
          }
          if(x_bound_2 == BoundaryCondition::Wall){
            // if(i1 == m_xdim-1 || i1 == m_xdim-2){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(i1 == xdim-1){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(x_bound_2 == BoundaryCondition::Open){
            if(i1 == xdim-1) i2 = i1;
          }
        }
        else{
          //Handle Y boundary conditions
          i0 = i1; i2 = i1;
          j0 = j1-1; j2 = j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+ydim)%ydim;
            j2 = (j2+ydim)%ydim;
          }
          if(y_bound_1 == BoundaryCondition::Wall){
            // if(j1 == 0 || j1 == 1){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(j1 == 0){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(y_bound_1 == BoundaryCondition::Open){
            if(j1 == 0) j0 = j1;
          }
          if(y_bound_2 == BoundaryCondition::Wall){
            // if(j1 == m_ydim-1 || j1 == m_ydim-2){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(j1 == ydim-1){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(y_bound_2 == BoundaryCondition::Open){
            if(j1 == ydim-1) j2 = j1;
          }
        }
        div(i1,j1) = (quantity(i2,j2) - quantity(i0,j0))/denom;
      }
    }
  }
  return div;
}

//Compute divergence term for simulation parameter "quantity"
//"quantity","vx","vy" used for transport term
//Non-transport terms contained in "nontransp_x", "nontransp_y"
Grid PlasmaDomain::divergence(const Grid &quantity, const Grid &nontransp_x, const Grid &nontransp_y, const Grid &vx, const Grid &vy){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result = transportDerivative1D(quantity, vx, 0) + transportDerivative1D(quantity, vy, 1);
  if(nontransp_x.size() > 1) result += derivative1D(nontransp_x, 0);
  if(nontransp_y.size() > 1) result += derivative1D(nontransp_y, 1);
  return result;
}


//Compute single-direction second derivative
Grid PlasmaDomain::secondDerivative1D(const Grid &quantity, const int index){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Zero(xdim,ydim);
  double denom = std::pow((m_dx*(1-index) + m_dy*(index)),2.0);
  #pragma omp parallel
  {
    #if BENCHMARKING_ON
    InstrumentationTimer timer("loop thread");
    #endif
    #pragma omp for collapse(2)
    for(int i=0; i<xdim; i++){
      for(int j=0; j<ydim; j++){
        int i0, i1, i2, j0, j1, j2;
        i1 = i; j1 = j;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j1; j2 = j1;
          i0 = i1-1; i2 = i1+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+xdim)%xdim;
            i2 = (i2+xdim)%xdim;
          }
          if(x_bound_1 == BoundaryCondition::Wall){
            // if(i1 == 0 || i1 == 1){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(i1 == 0){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(x_bound_1 == BoundaryCondition::Open){
            if(i1 == 0) i0 = i1;
          }
          if(x_bound_2 == BoundaryCondition::Wall){
            // if(i1 == m_xdim-1 || i1 == m_xdim-2){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(i1 == m_xdim-1){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(x_bound_2 == BoundaryCondition::Open){
            if(i1 == m_xdim-1) i2 = i1;
          }
        }
        else{
          //Handle Y boundary conditions
          i0 = i1; i2 = i1;
          j0 = j1-1; j2 = j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+ydim)%ydim;
            j2 = (j2+ydim)%ydim;
          }
          if(y_bound_1 == BoundaryCondition::Wall){
            // if(j1 == 0 || j1 == 1){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(j1 == 0){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(y_bound_1 == BoundaryCondition::Open){
            if(j1 == 0) j0 = j1;
          }
          if(y_bound_2 == BoundaryCondition::Wall){
            // if(j1 == m_ydim-1 || j1 == m_ydim-2){
            //   div(i1,j1) = 0.0;
            //   continue;
            // }
            if(j1 == m_ydim-1){
              div(i1,j1) = 0.0;
              continue;
            }
          }
          if(y_bound_2 == BoundaryCondition::Open){
            if(j1 == m_ydim-1) j2 = j1;
          }
        }
        div(i1,j1) = (quantity(i2,j2) - 2.0*quantity(i1,j1) + quantity(i0,j0))/denom;
      }
    }
  }
  return div;
}

//Computes Laplacian (del squared) of "quantity"
Grid PlasmaDomain::laplacian(const Grid &quantity){
  #if BENCHMARKING_ON
  InstrumentationTimer timer(__PRETTY_FUNCTION__);
  #endif
  Grid result_x = secondDerivative1D(quantity,0);
  Grid result_y = secondDerivative1D(quantity,1);
  return result_x+result_y;
}

