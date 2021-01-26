//mhdtoy
//Greg Szypko

#include <iostream>
#include <fstream>
#include <cmath>
#include "Eigen/Core"

#define XDIM 50
#define YDIM 50

#define DX 0.1
#define DY 0.1
#define NT 1000
#define OUTPUT_INTERVAL 5 //time steps between file outputs

#define EPSILON 0.1 //dynamic time stepping safety factor
#define GRIDFLOOR 0.0001 //min value for non-negative parameters
#define SUBCYCLES 10 //subcycles per cycle for thermal conduction

#define MU_0 1.0 //permeability of free space
#define K_B 1.0 //boltzmann constant
#define M_I 1.0 //ion mass
#define GRAV 1.0 //acceleration due to gravity
#define GAMMA 1.666667 //adiabatic index
#define KAPPA_0 1.0 //thermal conductivity coefficient
#define PI 3.14159265358979323846

//Boundary condition labels
#define PERIODIC 0
#define WALL 1
#define OPEN 2

#define XBOUND1 0
#define XBOUND2 0
#define YBOUND1 1
#define YBOUND2 2

//Output variables (0 to turn off output, 1 to turn on output)
#define RHO_OUT 1
#define TEMP_OUT 1
#define PRESS_OUT 0
#define ENERGY_OUT 0

using Grid = Eigen::Array<double,Eigen::Dynamic,Eigen::Dynamic>;

//Compute surface values from cell-centered values using Barton's method
//Meant to be used for transport terms only
//Result indexed s.t. element i,j indicates surface between i,j and i-1,j
//if "index"==0, or i,j and i,j-1 if "index"==1
Grid upwind_surface(const Grid &cell_center, const Grid &vel, const int index){
  int xdim = XDIM+1-index;
  int ydim = YDIM+index;
  Grid cell_surface = Grid::Zero(xdim,ydim);
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
        if(XBOUND1 == PERIODIC && XBOUND2 == PERIODIC){
          if(i2 == XDIM) continue;
          //Here, explicitly need macro XDIM instead of xdim
          i0 = (i0+XDIM)%XDIM;
          i1 = (i1+XDIM)%XDIM;
          i3 = (i3+XDIM)%XDIM;
        }
        if(XBOUND1 == WALL){
          if(i2 == 0){
            cell_surface(i2,j2) = 0.0;
            continue;
          } else if(i2 == 1){
            i0 = i1;
          }
        } else if(XBOUND1 == OPEN){
          if(i2 == 0){
            cell_surface(i2,j2) = 1.5*cell_center(i2,j2) - 0.5*cell_center(i3,j3); //lerp
            continue;
          } else if(i2 == 1){
            i0 = i1;
          }
        }
        if(XBOUND2 == WALL){
          if(i2 == XDIM){
            cell_surface(i2,j2) = 0.0;
            continue;
          } else if(i2 == XDIM - 1){
            i3 = i2;
          }
        } else if(XBOUND2 == OPEN){
          if(i2 == XDIM){
            cell_surface(i2,j2) = 1.5*cell_center(i1,j1) - 0.5*cell_center(i0,j0); //lerp
            continue;
          } else if(i2 == XDIM - 1){
            i3 = i2;
          }
        }
      }
      else{
        //Handle Y boundary conditions
        i0 = i2; i1 = i2; i3 = i2;
        j0 = j2-2; j1 = j2-1; j3 = j2+1;
        if(YBOUND1 == PERIODIC && YBOUND2 == PERIODIC){
          if(j2 == YDIM) continue;
          //Here, explicitly need macro YDIM instead of ydim
          j0 = (j0+YDIM)%YDIM;
          j1 = (j1+YDIM)%YDIM;
          j3 = (j3+YDIM)%YDIM;
        }
        if(YBOUND1 == WALL){
          if(j2 == 0){
            cell_surface(i2,j2) = 0.0;
            continue;
          } else if(j2 == 1){
            j0 = j1;
          }
        } else if(YBOUND1 == OPEN){
          if(j2 == 0){
            cell_surface(i2,j2) = 1.5*cell_center(i2,j2) - 0.5*cell_center(i3,j3); //lerp
            continue;
          } else if(j2 == 1){
            j0 = j1;
          }
        }
        if(YBOUND2 == WALL){
          if(j2 == YDIM){
            cell_surface(i2,j2) = 0.0;
            continue;
          } else if(j2 == YDIM - 1){
            j3 = j2;
          }
        } else if(YBOUND2 == OPEN){
          if(j2 == YDIM){
            cell_surface(i2,j2) = 1.5*cell_center(i1,j1) - 0.5*cell_center(i0,j0); //lerp
            continue;
          } else if(j2 == YDIM - 1){
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
  return cell_surface;
}

Grid transport_divergence1D(const Grid &quantity, const Grid &vel, const int index){
  Grid surf_quantity = upwind_surface(quantity, vel, index);
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  double denom = DX*(1-index) + DY*(index);
  Grid div = Grid::Ones(xdim,ydim);
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      int i0, i1, i2, i2surf, j0, j1, j2, j2surf; //Need separate indices for surface and vel
      i1 = i; j1 = j;
      if(index == 0){
        //Handle X boundary conditions
        j0 = j1; j2 = j1; j2surf = j1;
        i0 = i1-1; i2 = i1+1; i2surf = i1+1;
        if(XBOUND1 == PERIODIC && XBOUND2 == PERIODIC){
          i0 = (i0+XDIM)%XDIM;
          i2 = (i2+XDIM)%XDIM;
          i2surf = (i2surf+XDIM)%XDIM;
        }
        if(XBOUND1 == WALL || XBOUND1 == OPEN){
          if(i1 == 0) i0 = i1;
        }
        if(XBOUND2 == WALL || XBOUND2 == OPEN){
          if(i1 == XDIM-1) i2 = i1;
        }
      }
      else{
        //Handle Y boundary conditions
        i0 = i1; i2 = i1; i2surf = i1;
        j0 = j1-1; j2 = j1+1; j2surf = j1+1;
        if(YBOUND1 == PERIODIC && YBOUND2 == PERIODIC){
          j0 = (j0+YDIM)%YDIM;
          j2 = (j2+YDIM)%YDIM;
          j2surf = (j2surf+YDIM)%YDIM;
        }
        if(YBOUND1 == WALL || YBOUND1 == OPEN){
          if(j1 == 0) j0 = j1;
        }
        if(YBOUND2 == WALL || YBOUND2 == OPEN){
          if(j1 == YDIM-1) j2 = j1;
        }
      }
      div(i1,j1) = (surf_quantity(i2surf,j2surf)*0.5*(vel(i1,j1)+vel(i2,j2))
                - surf_quantity(i1,j1)*0.5*(vel(i1,j1)+vel(i0,j0)))/denom;
    }
  }
  return div;
}


//Compute single-direction divergence term for non-transport term (central differencing)
Grid divergence1D(const Grid &quantity, const int index){
  int xdim = quantity.rows();
  int ydim = quantity.cols();
  Grid div = Grid::Ones(xdim,ydim);
  double denom = 2.0*(DX*(1-index) + DY*(index));
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      int i0, i1, i2, j0, j1, j2;
      i1 = i; j1 = j;
      if(index == 0){
        //Handle X boundary conditions
        j0 = j1; j2 = j1;
        i0 = i1-1; i2 = i1+1;
        //ENFORCES PERIODIC X-BOUNDARIES
        if(XBOUND1 == PERIODIC && XBOUND2 == PERIODIC){
          i0 = (i0+xdim)%xdim;
          i2 = (i2+xdim)%xdim;
        }
        if(XBOUND1 == WALL || XBOUND1 == OPEN){
          if(i1 == 0) i0 = i1;
        }
        if(XBOUND2 == WALL || XBOUND2 == OPEN){
          if(i1 == XDIM-1) i2 = i1;
        }
      }
      else{
        //Handle Y boundary conditions
        i0 = i1; i2 = i1;
        j0 = j1-1; j2 = j1+1;
        if(YBOUND1 == PERIODIC && YBOUND2 == PERIODIC){
          j0 = (j0+ydim)%ydim;
          j2 = (j2+ydim)%ydim;
        }
        if(YBOUND1 == WALL || YBOUND1 == OPEN){
          if(j1 == 0) j0 = j1;
        }
        if(YBOUND2 == WALL || YBOUND2 == OPEN){
          if(j1 == YDIM-1) j2 = j1;
        }
      }
      div(i1,j1) = (quantity(i2,j2) - quantity(i0,j0))/denom;
    }
  }
  return div;
}

//Compute divergence term for simulation parameter "quantity"
//"quantity","vx","vy" used for transport term
//Non-transport terms contained in "nontransp_x", "nontransp_y"
Grid divergence(const Grid &quantity, const Grid &nontransp_x, const Grid &nontransp_y, const Grid &vx, const Grid &vy){
  // Eigen::IOFormat one_line_format(4, Eigen::DontAlignCols, ",", ";", "", "", "", "\n");
  Grid result = transport_divergence1D(quantity, vx, 0) + transport_divergence1D(quantity, vy, 1);
  if(nontransp_x.size() > 1) result += divergence1D(nontransp_x, 0);
  if(nontransp_y.size() > 1) result += divergence1D(nontransp_y, 1);
  // std::cout << result.matrix().format(one_line_format);
  return result;
}

//Enforce dynamic time stepping
double recompute_dt(const Grid &press, const Grid &rho, const Grid &vx, const Grid &vy){
  double running_min_dt = std::numeric_limits<double>::max();
  for(int i=0; i<XDIM; i++){
    for(int j=0; j<YDIM; j++){
      double c_s = std::sqrt(GAMMA*press(i,j)/rho(i,j));
      double abs_vx = DX/(c_s+std::abs(vx(i,j)));
      double abs_vy = DY/(c_s+std::abs(vy(i,j)));
      running_min_dt = std::min(running_min_dt, abs_vx);
      running_min_dt = std::min(running_min_dt, abs_vy);
    }
  }
  return EPSILON*running_min_dt;
}

//Generates gaussian initial condition for a variable, centered at middle of grid
Grid GaussianGrid(int xdim, int ydim, double min, double max){
  Eigen::VectorXd gauss_x(xdim), gauss_y(ydim);
  double sigmax = 0.05*xdim;
  double sigmay = 0.05*ydim;
  for(int i=0; i<xdim; i++){
    gauss_x(i) = std::exp(-0.5*std::pow(((double)i-0.5*(double)(xdim-1))/sigmax,2.0));
  }
  for(int j=0; j<ydim; j++){
    gauss_y(j) = std::exp(-0.5*std::pow(((double)j-0.5*(double)(ydim-1))/sigmay,2.0));
  }
  Eigen::MatrixXd gauss_2d = gauss_x * gauss_y.transpose();
  gauss_2d = gauss_2d*(max-min) + min*Eigen::MatrixXd::Ones(xdim,ydim);
  return gauss_2d.array();
}

//Generates potential bipolar field for component corresponding to index "index"
//Centered s.t. origin lies at bottom middle of domain
//Pressure scale height h, field poles at +/- l, field strength at poles b0
Grid BipolarField(const int xdim, const int ydim, const double b0, const double h, const int index){
  Grid result = Grid::Zero(xdim, ydim);
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      double x = (i - (double)(xdim-1)*0.5)*DX;
      double y = (j - (double)(ydim-1)*0.5)*DY;
      if(index == 0) result(i,j) = std::exp(-0.5*y/h)*std::cos(0.5*x/h);
      else result(i,j) = -std::exp(-0.5*y/h)*std::sin(0.5*x/h);
    }
  }
  return result;
}

//Computes cell-centered conductive flux from temperature "temp"
//Flux computed in direction indicated by "index": 0 for x, 1 for y
//k0 is conductive coefficient
Grid conductive_flux(const Grid &temp, const double k0, const int index){
  int xdim = temp.rows();
  int ydim = temp.cols();
  Grid flux = Grid::Zero(xdim,ydim);
  for(int i=0; i<xdim; i++){
    for(int j=0; j<ydim; j++){
      flux(i,j) = std::pow(temp(i,j),7.0/2.0);
    }
  }
  return -2.0/7.0*k0*divergence1D(flux,index);
}

int main(int argc,char* argv[]){
  Eigen::IOFormat one_line_format(4, Eigen::DontAlignCols, ",", ";", "", "", "", "\n");
  std::ofstream out_file;
  out_file.open ("output.txt");
  out_file << XDIM << "," << YDIM << std::endl;

  Grid rho, mom_x, mom_y, temp, press, energy, bx, by, bz;

  // rho = Grid::Ones(XDIM,YDIM);
  // for(int i=0; i<XDIM; i++)
  // for(int j=0; j<YDIM; j++)
  // rho(i,j) = (1.0 - 0.1*abs(i-XDIM*0.5)/XDIM);

  rho = GaussianGrid(XDIM, YDIM, 1.0, 5.0);

  mom_x = Grid::Zero(XDIM,YDIM); //x momentum density
  mom_y = Grid::Zero(XDIM,YDIM); //y momentum density
  // temp = Grid::Ones(XDIM,YDIM); //temperature
  temp = GaussianGrid(XDIM, YDIM, 1.0, 2.0); //temperature
  double b0 = 1.0;
  double h = 4.0*PI/(XDIM*DX);
  bx = BipolarField(XDIM, YDIM, b0, h, 0);
  by = BipolarField(XDIM, YDIM, b0, h, 1);
  bz = Grid::Zero(XDIM,YDIM);
  // bx = Grid::Zero(XDIM,YDIM);
  // by = Grid::Zero(XDIM,YDIM);
  // bz = Grid::Zero(XDIM,YDIM);

  out_file << bx.matrix().format(one_line_format);
  out_file << by.matrix().format(one_line_format);

  Grid mag_press = (bx*bx + by*by + bz*bz)/(2.0*MU_0);
  Grid mag_pxx = (-bx*bx + by*by + bz*bz)/(2.0*MU_0);
  Grid mag_pyy = (bx*bx - by*by + bz*bz)/(2.0*MU_0);
  Grid mag_pzz = (bx*bx + by*by - bz*bz)/(2.0*MU_0);
  Grid mag_pxy = -bx*by/MU_0;
  Grid mag_pxz = -bx*bz/MU_0;
  Grid mag_pyz = -by*bz/MU_0;
  press = 2.0*K_B*rho*temp/M_I; //assumes 
  energy = press/(GAMMA - 1.0) + 0.5*(mom_x*mom_x + mom_y*mom_y)/rho + mag_press;

  //Output intial state
  double t=0.0;
  double dt=1.0;
  out_file << "t=" << t << std::endl;
  if(RHO_OUT){
    out_file << "rho\n";
    out_file << rho.matrix().format(one_line_format);
  }
  if(TEMP_OUT){
    out_file << "temp\n";
    out_file << temp.matrix().format(one_line_format);
  }
  if(PRESS_OUT){
    out_file << "press\n";
    out_file << press.matrix().format(one_line_format);
  }
  if(ENERGY_OUT){
    out_file << "energy\n";
    out_file << energy.matrix().format(one_line_format);
  }

  for (int iter = 0; iter < NT; iter++){
    // Enforce rigid lower boundary
    // for(int i=0; i<XDIM; i++){
    //   mom_x(i,0) = 0.0;
    //   mom_y(i,0) = 0.0;
    // }

    //Compute values needed for time evolution
    Grid vx = mom_x/rho;
    Grid vy = mom_y/rho;
    press = 2.0*K_B*rho*temp/M_I; //assumes 
    energy = press/(GAMMA - 1.0) + 0.5*(mom_x*vx + mom_y*vy) + mag_press;
    // Grid press = (GAMMA - 1.0)*(energy - 0.5*(mom_x*vx + mom_y*vy));
    dt = recompute_dt(press, rho, vx, vy);

    //Subcycle to simulate thermal diffusion
    Grid energy_relaxed = energy;
    for(int subcycle = 0; subcycle < SUBCYCLES; subcycle++){
      Grid con_flux_x = conductive_flux(temp, KAPPA_0, 0);
      Grid con_flux_y = conductive_flux(temp, KAPPA_0, 1);
      energy_relaxed = energy_relaxed - (dt/(double)SUBCYCLES)*(divergence1D(con_flux_x,0)+divergence1D(con_flux_y,1));
      press = (GAMMA - 1.0)*(energy_relaxed - 0.5*(mom_x*vx + mom_y*vy) - mag_press);
      temp = M_I*press/(2.0*K_B*rho);
    }

    // std::cout << "rho: " << rho.matrix().format(one_line_format);
    // std::cout << "mom_x: " << mom_x.matrix().format(one_line_format);
    // std::cout << "mom_y: " << mom_y.matrix().format(one_line_format);
    // std::cout << "energy: " << energy.matrix().format(one_line_format);
    // std::cout << "press: " << press.matrix().format(one_line_format);
    //Advance time by dt
    Grid zero = Grid::Zero(1,1);
    Grid rho_next = rho - dt*divergence(rho,zero,zero,vx,vy);
    Grid mom_x_next = mom_x - dt*divergence(mom_x, press + mag_pxx, mag_pxy, vx, vy);
    Grid mom_y_next = mom_y - dt*divergence(mom_y, mag_pxy, press + mag_pyy, vx, vy) - dt*rho*GRAV;
    //NOTE: This is a naive implementation of thermal conduction; need to implement subcycling for thermal diffusion
    Grid energy_next = energy_relaxed - dt*divergence(energy_relaxed+press, mag_pxx*vx + mag_pxy*vy, mag_pxy*vx + mag_pyy*vy, vx, vy) - dt*rho*vy*GRAV;
    //CENTRAL DIFFERENCING VERSION
    // Grid rho_next = rho - dt*divergence(zero,rho*vx,rho*vy,vx,vy);
    // Grid mom_x_next = mom_x - dt*divergence(zero, mom_x*vx+press, mom_x*vy, vx, vy);
    // Grid mom_y_next = mom_y - dt*divergence(zero, mom_y*vx, mom_y*vy+press, vx, vy);
    // Grid energy_next = energy - dt*divergence(zero, energy*vx+press*vx, energy*vy+press*vy, vx, vy);

    //Sanity checks
    for(int i=0; i<XDIM; i++){
      for(int j=0; j<YDIM; j++){
        //Ensure no negative densities
        rho_next(i,j) = std::max(GRIDFLOOR, rho_next(i,j));
        //Ensure no negative pressure/temperature
        double bulk_kin_energy = 0.5*(mom_x_next(i,j)*mom_x_next(i,j)
                                    + mom_y_next(i,j)*mom_y_next(i,j))/rho_next(i,j);
        energy_next(i,j) = std::max(bulk_kin_energy + mag_press(i,j), energy_next(i,j));
      }
    }
    //Swap current and next pointers
    rho = rho_next;
    mom_x = mom_x_next;
    mom_y = mom_y_next;
    press = (GAMMA - 1.0)*(energy_next - 0.5*(mom_x*vx + mom_y*vy) - mag_press);
    temp = M_I*press/(2.0*K_B*rho);

    t = t + dt;
    if(iter%10 == 0){
      printf("\rIterations: %i/%i", iter, NT);
      fflush(stdout);
    }

    if(iter%OUTPUT_INTERVAL == 0){
      out_file << "t=" << t << std::endl;
      if(RHO_OUT){
        out_file << "rho\n";
        out_file << rho.matrix().format(one_line_format);
      }
      if(TEMP_OUT){
        out_file << "temp\n";
        out_file << temp.matrix().format(one_line_format);
      }
      if(PRESS_OUT){
        out_file << "press\n";
        out_file << press.matrix().format(one_line_format);
      }
      if(ENERGY_OUT){
        out_file << "energy\n";
        out_file << energy.matrix().format(one_line_format);
      }
    }
  }
  std::cout << "\rIterations: " << NT << "/" << NT << "\n";
  out_file.close();
}