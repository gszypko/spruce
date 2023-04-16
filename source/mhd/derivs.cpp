//derivs.cpp
//Defines differentiation functions for PlasmaDomain

#include "plasmadomain.hpp"

//Compute surface values from cell-centered values using Barton's method
//Meant to be used for transport terms only
//Result indexed s.t. element i,j indicates surface between i,j and i-1,j
//if "index"==0, or i,j and i,j-1 if "index"==1
Grid PlasmaDomain::upwindSurface(const Grid &cell_center, const Grid &vel, const int index, int xl, int yl, int xu, int yu){
  // if(index == 0) assert(xl>1 && xu<cell_center.rows()-2 && "Must have at least two-cell border in direction of differentiation for Barton's method");
  // else if(index == 1) assert(yl>1 && yu<cell_center.cols()-2 && "Must have at least two-cell border in direction of differentiation for Barton's method");
  int xdim = m_xdim+1-index;
  int ydim = m_ydim+index;
  Grid cell_surface = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu+1-index; i++){
      for(int j = yl; j <= yu+index; j++){
        //Handle direction of cell_center being considered (i.e. index for differencing)
        int i2 = i, j2 = j;
        int i0, i1, i3, j0, j1, j3;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j2; j1 = j2; j3 = j2;
          i0 = i2-2; i1 = i2-1; i3 = i2+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            if(i2 == m_xdim) continue; //So we don't double count the periodic boundary face
            i0 = (i0+m_xdim)%m_xdim;
            i1 = (i1+m_xdim)%m_xdim;
            i3 = (i3+m_xdim)%m_xdim;
          }
        }
        else{
          //Handle Y boundary conditions
          i0 = i2; i1 = i2; i3 = i2;
          j0 = j2-2; j1 = j2-1; j3 = j2+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            if(j2 == m_ydim) continue; //So we don't double count the periodic boundary face
            j0 = (j0+m_ydim)%m_ydim;
            j1 = (j1+m_ydim)%m_ydim;
            j3 = (j3+m_ydim)%m_ydim;
          }
        }
        //Apply Barton's method
        double d1, d2, d3;
        d2 = boundaryInterpolate(cell_center,i1,j1,i2,j2);
        if(boundaryInterpolate(vel,i1,j1,i2,j2)>0.0){
          d3 = cell_center(i1,j1);
          d1 = boundaryExtrapolate(cell_center,i0,j0,i1,j1);
          if(cell_center(i2,j2) <= cell_center(i1,j1)){
            cell_surface(i2,j2) = std::min(d3,std::max(d1,d2));
          } else { //cell_center(i2,j2) > cell_center(i1,j1)
            cell_surface(i2,j2) = std::max(d3,std::min(d1,d2));
          }
        } else if(boundaryInterpolate(vel,i1,j1,i2,j2)<0.0){
          d3 = cell_center(i2,j2);
          d1 = boundaryExtrapolate(cell_center,i3,j3,i2,j2);
          if(cell_center(i2,j2) <= cell_center(i1,j1)){
            cell_surface(i2,j2) = std::max(d3,std::min(d1,d2));
          } else { //cell_center(i2,j2) > cell_center(i1,j1)
            cell_surface(i2,j2) = std::min(d3,std::max(d1,d2));
          }
        } else {
          cell_surface(i2,j2) = d2; //default case if no flow direction; linear interp
        }
      }
    }
  }
  return cell_surface;
}

// Specialized version of PlasmaDomain::upwindSurface, used when the advection direction is uniform and
// known ahead of time. Requires only three-cell template instead of the four-cell template for upwindSurface (which permits bidirectional velocity).
Grid PlasmaDomain::upwindSurfaceDirected(const Grid &cell_center, const bool vel_positive, const int index, int xl, int yl, int xu, int yu){
  int xdim = m_xdim+1-index;
  int ydim = m_ydim+index;
  Grid cell_surface = Grid::Zero(xdim,ydim);
  int vel_sign = (vel_positive ? 1 : -1);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
        //Handle direction of cell_center being considered (i.e. index for differencing)
        int i2 = i, j2 = j;
        int i0, i1, j0, j1;
        if(index == 0){
          j0 = j2; j1 = j2;
          i0 = i2-(vel_sign*2); i1 = i2-(vel_sign*1);
          assert(x_bound_1 != BoundaryCondition::Periodic && x_bound_2 != BoundaryCondition::Periodic && "This should never be used with periodic boundaries");
        }
        else{
          i0 = i2; i1 = i2;
          j0 = j2-(vel_sign*2); j1 = j2-(vel_sign*1);
          assert(y_bound_1 != BoundaryCondition::Periodic && y_bound_2 != BoundaryCondition::Periodic && "This should never be used with periodic boundaries");
        }
        //Apply Barton's method
        double d1, d2, d3;
        d2 = boundaryInterpolate(cell_center,i1,j1,i2,j2);
        d3 = cell_center(i1,j1);
        d1 = boundaryExtrapolate(cell_center,i0,j0,i1,j1);
        int i_surface, j_surface; //surface quantities returned s.t. (i,j) corresponds to the surface between cells (i-1,j) and (i,j) or (i,j-1) and (i,j), depending on the index
        i_surface = i2 + (index == 0 && !vel_positive ? 1 : 0);
        j_surface = j2 + (index == 1 && !vel_positive ? 1 : 0);
        if(cell_center(i2,j2) <= cell_center(i1,j1)){
          cell_surface(i_surface,j_surface) = std::min(d3,std::max(d1,d2));
        } else { //cell_center(i2,j2) > cell_center(i1,j1)
          cell_surface(i_surface,j_surface) = std::max(d3,std::min(d1,d2));
        }
      }
    }
  }
  return cell_surface;
}

// Compute 1D derivative of quantity using Barton's method for upwinding (based on the velocity vel)
// When multiply_vel is true, computes derivative of (quantity*vel) -- to be used for advection term of conserved quantity, i.e. a flux divergence
// When multiply_vel is false, computes derivative of quantity while using vel only for upwinding direction
Grid PlasmaDomain::transportDerivative1D(const Grid &quantity, const Grid &vel, const int index, int xl, int yl, int xu, int yu, bool multiply_vel){
  // if(index == 0) assert(xl>0 && xu<quantity.rows()-1 && "Must have at least one-cell border in direction of differentiation");
  // else if(index == 1) assert(yl>0 && yu<quantity.cols()-1 && "Must have at least one-cell border in direction of differentiation");
  Grid surf_quantity = upwindSurface(quantity, vel, index, xl, yl, xu, yu);
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid denom = m_grids[d_x]*(double)(1-index) + m_grids[d_y]*(double)(index);
  Grid div = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
        int i0, i1, i2, j0, j1, j2;
        i1 = i; j1 = j;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j1; j2 = j1; //j2surf = j1;
          i0 = i1-1; i2 = i1+1; //i2surf = i1+1;
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+m_xdim)%m_xdim;
            i2 = (i2+m_xdim)%m_xdim;
          }
        }
        else{
          //Handle Y boundary conditions
          i0 = i1; i2 = i1; //i2surf = i1;
          j0 = j1-1; j2 = j1+1; //j2surf = j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+m_ydim)%m_ydim;
            j2 = (j2+m_ydim)%m_ydim;
          }
        }
        if(multiply_vel) div(i1,j1) = (surf_quantity(i2,j2)*boundaryInterpolate(vel,i1,j1,i2,j2)
                  - surf_quantity(i1,j1)*boundaryInterpolate(vel,i0,j0,i1,j1))/denom(i1,j1);
        else div(i1,j1) = (surf_quantity(i2,j2) - surf_quantity(i1,j1))/denom(i1,j1);
      }
    }
  }
  return div;
}

// Similar to transportDerivative1D, but does not take a velocity term for upwind determination
// This function meant to be used when the upwind direction can be determined/enforced ahead of time,
// as is the case for boundary-perpendicular gradients (that aren't negated) in the method of characteristics
Grid PlasmaDomain::characteristicBartonDerivative1D(const Grid &quantity, const bool positive_forward, const int index, int xl, int yl, int xu, int yu){
  Grid surf_quantity = upwindSurfaceDirected(quantity, positive_forward, index, xl, yl, xu, yu);
  // Fill in outer edge surface quantites (as direct upwind quantity)
  if(index == 0){
    if(positive_forward){
      for(int j=yl; j<=yu; j++) surf_quantity(surf_quantity.rows()-1,j) = quantity(quantity.rows()-1,j)
                                                                        + (quantity(quantity.rows()-1,j) - surf_quantity(surf_quantity.rows()-2,j));
    } else{
      for(int j=yl; j<=yu; j++) surf_quantity(0,j) = quantity(0,j)
                                                      - (surf_quantity(1,j) - quantity(0,j));
    }
  }else {
    assert(index == 1);
    if(positive_forward){
      for(int i=xl; i<=xu; i++) surf_quantity(i,surf_quantity.cols()-1) = quantity(i,quantity.cols()-1)
                                                                          + (quantity(i,quantity.cols()-1) - surf_quantity(i,surf_quantity.cols()-2));
    } else{
      for(int i=xl; i<=xu; i++) surf_quantity(i,0) = quantity(i,0)
                                                    - (surf_quantity(i,1) - quantity(i,0));
    }
  }
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid denom = m_grids[d_x]*(double)(1-index) + m_grids[d_y]*(double)(index);
  Grid div = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
        int i1, i2, j1, j2;
        i1 = i; j1 = j;
        if(index == 0){
          j2 = j1;
          i2 = i1+1;
          assert(x_bound_1 != BoundaryCondition::Periodic && x_bound_2 != BoundaryCondition::Periodic && "This should never be used with periodic boundaries");
        }
        else{
          i2 = i1;
          j2 = j1+1;
          assert(y_bound_1 != BoundaryCondition::Periodic && y_bound_2 != BoundaryCondition::Periodic && "This should never be used with periodic boundaries");
        }
        div(i1,j1) = (surf_quantity(i2,j2) - surf_quantity(i1,j1))/denom(i1,j1);
      }
    }
  }
  return div;
}

Grid PlasmaDomain::transportDivergence2D(const Grid &quantity, const std::vector<Grid> &vel, int xl, int yl, int xu, int yu)
{
  assert(vel.size() == 2 && "this operator requires 2D velocity");
  return transportDerivative1D(quantity,vel[0],0, xl, yl, xu, yu) + transportDerivative1D(quantity,vel[1],1, xl, yl, xu, yu);
}

//Compute single-direction divergence term for non-transport term (central differencing)
Grid PlasmaDomain::derivative1D(const Grid &quantity, const int index, int xl, int yl, int xu, int yu){
  // if(index == 0) assert(xl>0 && xu<quantity.rows()-1 && "Must have at least one-cell border in direction of differentiation");
  // else if(index == 1) assert(yl>0 && yu<quantity.cols()-1 && "Must have at least one-cell border in direction of differentiation");
  // else assert(false && "This function assumes two dimensions");
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
        int i0, i1, i2, j0, j1, j2;
        i1 = i; j1 = j;
        double denom;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j1; j2 = j1;
          i0 = i1-1; i2 = i1+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+xdim)%xdim;
            i2 = (i2+xdim)%xdim;
          }
          denom = m_grids[d_x](i1,j1);
        }
        else{
          //Handle Y boundary conditions
          i0 = i1; i2 = i1;
          j0 = j1-1; j2 = j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+ydim)%ydim;
            j2 = (j2+ydim)%ydim;
          }
          denom = m_grids[d_y](i1,j1);
        }
        div(i1,j1) = (boundaryInterpolate(quantity,i1,j1,i2,j2) - boundaryInterpolate(quantity,i0,j0,i1,j1)) / denom;
      }
    }
  }
  return div;
}

Grid PlasmaDomain::derivative1DBackward2ndOrder(const Grid &quantity, bool positive_forward, const int index, int xl, int yl, int xu, int yu){
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
        int i0, i1, i2, j0, j1, j2;
        i2 = i; j2 = j;
        // double denom;
        double delta1, delta2, deltafull;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j1 = j2;
          i1 = positive_forward ? i2-1 : i2+1;
          i0 = positive_forward ? i1-1 : i1+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+xdim)%xdim;
            i1 = (i1+xdim)%xdim;
          }
          delta1 = 0.5*(m_grids[d_x](i0,j0) + m_grids[d_x](i1,j1));
          delta2 = 0.5*(m_grids[d_x](i1,j1) + m_grids[d_x](i2,j2));
          // denom = 0.5*(m_grids[d_x](i1,j1) + m_grids[d_x](i0,j0));
        }
        else{
          //Handle Y boundary conditions
          i0 = i1 = i2;
          j1 = positive_forward ? j2-1 : j2+1;
          j0 = positive_forward ? j1-1 : j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+ydim)%ydim;
            j1 = (j1+ydim)%ydim;
          }
          delta1 = 0.5*(m_grids[d_y](i0,j0) + m_grids[d_y](i1,j1));
          delta2 = 0.5*(m_grids[d_y](i1,j1) + m_grids[d_y](i2,j2));
          // denom = 0.5*(m_grids[d_y](i1,j1) + m_grids[d_y](i0,j0));
        }
        // div(i1,j1) = (positive_forward ? 1.0 : -1.0)*(quantity(i1,j1) - quantity(i0,j0)) / denom;
        deltafull = delta1 + delta2;
        div(i2,j2) = (positive_forward ? 1.0 : -1.0)*(
          (deltafull + delta2)/(delta2*deltafull)*quantity(i2,j2) - deltafull/(delta1*delta2)*quantity(i1,j1) + delta2/(deltafull*delta1)*quantity(i0,j0)
        );
      }
    }
  }
  
  return div;
}

// First-order backward difference
Grid PlasmaDomain::derivative1DBackward(const Grid &quantity, bool positive_forward, const int index, int xl, int yl, int xu, int yu){
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
        int i0, i1, j0, j1;
        i1 = i; j1 = j;
        double denom;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j1;
          i0 = positive_forward ? i1-1 : i1+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+xdim)%xdim;
          }
          denom = 0.5*(m_grids[d_x](i1,j1) + m_grids[d_x](i0,j0));
        }
        else{
          //Handle Y boundary conditions
          i0 = i1;
          j0 = positive_forward ? j1-1 : j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+ydim)%ydim;
          }
          denom = 0.5*(m_grids[d_y](i1,j1) + m_grids[d_y](i0,j0));
        }
        div(i1,j1) = (positive_forward ? 1.0 : -1.0)*(quantity(i1,j1) - quantity(i0,j0)) / denom;
      }
    }
  }
  
  return div;
}

Grid PlasmaDomain::secondDerivative1DBackward(const Grid &quantity, bool positive_forward, const int index, int xl, int yl, int xu, int yu){
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Zero(xdim,ydim);
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
        int i0, i1, i2, j0, j1, j2;
        i2 = i; j2 = j;
        double denom01, denom12;
        if(index == 0){
          //Handle X boundary conditions
          j0 = j2; j1 = j2;
          i1 = positive_forward ? i2-1 : i2+1;
          i0 = positive_forward ? i1-1 : i1+1;
          //ENFORCES PERIODIC X-BOUNDARIES
          if(x_bound_1 == BoundaryCondition::Periodic && x_bound_2 == BoundaryCondition::Periodic){
            i0 = (i0+xdim)%xdim;
            i1 = (i1+xdim)%xdim;
          }
          denom01 = 0.5*(m_grids[d_x](i1,j1) + m_grids[d_x](i0,j0));
          denom12 = 0.5*(m_grids[d_x](i1,j1) + m_grids[d_x](i2,j2));
        }
        else{
          //Handle Y boundary conditions
          i0 = i2; i1 = i2;
          j1 = positive_forward ? j2-1 : j2+1;
          j0 = positive_forward ? j1-1 : j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+ydim)%ydim;
            j1 = (j1+ydim)%ydim;
          }
          denom01 = 0.5*(m_grids[d_y](i1,j1) + m_grids[d_y](i0,j0));
          denom12 = 0.5*(m_grids[d_y](i1,j1) + m_grids[d_y](i2,j2));
        }
        // div(i1,j1) = (quantity(std::max(i1,i0),std::max(j1,j0)) - quantity(std::min(i1,i0),std::min(j1,j0))) / denom;
        div(i2,j2) = (positive_forward ? 1.0 : -1.0)*(
                        (quantity(std::max(i1,i2),std::max(j1,j2)) - quantity(std::min(i1,i2),std::min(j1,j2))) / denom12
                      - (quantity(std::max(i1,i0),std::max(j1,j0)) - quantity(std::min(i1,i0),std::min(j1,j0))) / denom01) / denom12;
      }
    }
  }
  
  return div;
}


Grid PlasmaDomain::divergence2D(const Grid& a_x, const Grid& a_y, int xl, int yl, int xu, int yu){
  return derivative1D(a_x, 0, xl, yl, xu, yu) + derivative1D(a_y, 1, xl, yl, xu, yu);
}

Grid PlasmaDomain::divergence2D(const std::vector<Grid>& a, int xl, int yl, int xu, int yu){
  assert(a.size() == 2 && "divergence function assumes two vector components");
  return divergence2D(a[0],a[1], xl, yl, xu, yu);
}

//Compute single-direction second derivative
Grid PlasmaDomain::secondDerivative1D(const Grid &quantity, const int index, int xl, int yl, int xu, int yu){
  // if(index == 0) assert(xl>0 && xu<quantity.rows()-1 && "Must have at least one-cell border in direction of differentiation");
  // else if(index == 1) assert(yl>0 && yu<quantity.cols()-1 && "Must have at least one-cell border in direction of differentiation");
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Zero(xdim,ydim);
  Grid denom = (0.5*m_grids[d_x]*(double)(1-index) + 0.5*m_grids[d_y]*(double)(index)).square();
  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = xl; i <= xu; i++){
      for(int j = yl; j <= yu; j++){
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
        }
        else{
          //Handle Y boundary conditions
          i0 = i1; i2 = i1;
          j0 = j1-1; j2 = j1+1;
          if(y_bound_1 == BoundaryCondition::Periodic && y_bound_2 == BoundaryCondition::Periodic){
            j0 = (j0+ydim)%ydim;
            j2 = (j2+ydim)%ydim;
          }
        }
        div(i1,j1) = (boundaryInterpolate(quantity,i1,j1,i2,j2) - 2.0*quantity(i1,j1) + boundaryInterpolate(quantity,i0,j0,i1,j1))/denom(i1,j1);
      }
    }
  }
  return div;
}

//Computes Laplacian (del squared) of "quantity"
Grid PlasmaDomain::laplacian(const Grid &quantity, int xl, int yl, int xu, int yu){
  Grid result_x = secondDerivative1D(quantity,0,xl,yl,xu,yu);
  Grid result_y = secondDerivative1D(quantity,1,xl,yl,xu,yu);
  return result_x+result_y;
}

//Computes curl of vector in z-direction (result in xy-plane)
std::vector<Grid> PlasmaDomain::curlZ(const Grid& z, int xl, int yl, int xu, int yu){
  Grid result_x = derivative1D(z,1,xl,yl,xu,yu);
  Grid result_y = -derivative1D(z,0,xl,yl,xu,yu);
  return {result_x, result_y};
}

//Computes curl of vector in xy-plane (result in z-direction)
Grid PlasmaDomain::curl2D(const Grid& x, const Grid& y, int xl, int yl, int xu, int yu){
  return derivative1D(y, 0,xl,yl,xu,yu) - derivative1D(x, 1,xl,yl,xu,yu);
}

//Linearly interpolate value of quantity for boundary between (i1,j1) and (i2,j2)
double PlasmaDomain::boundaryInterpolate(const Grid &quantity, int i1, int j1, int i2, int j2){
  assert(i1 == i2 || j1 == j2 && "boundary interpolation only meant to operate along one axis");
  double a = quantity(i1,j1), b = quantity(i2,j2);
  double dist_a, dist_b;
  if(i1 == i2){
    dist_a = 0.5*m_grids[d_y](i1,j1); dist_b = 0.5*m_grids[d_y](i2,j2);
  } else { //j1 == j2
    dist_a = 0.5*m_grids[d_x](i1,j1); dist_b = 0.5*m_grids[d_x](i2,j2);
  }
  return (a*dist_b + b*dist_a)/(dist_b + dist_a);
}

//Linearly extrapolate value of quantity for between (i2,j2) and (i3,j3) using (i1,j1) and (i2,j2)
double PlasmaDomain::boundaryExtrapolate(const Grid &quantity, int i1, int j1, int i2, int j2){
  assert(i1 == i2 || j1 == j2 && "boundary interpolation only meant to operate along one axis");
  double a = quantity(i1,j1), b = quantity(i2,j2);
  double dist_a, dist_b;
  if(i1 == i2){
    dist_a = 0.5*m_grids[d_y](i1,j1); dist_b = 0.5*m_grids[d_y](i2,j2);
  } else { //j1 == j2
    dist_a = 0.5*m_grids[d_x](i1,j1); dist_b = 0.5*m_grids[d_x](i2,j2);
  }
  return a + (b-a)*(dist_a + 2.0*dist_b)/(dist_a + dist_b);
}