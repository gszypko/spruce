#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

/***************************** SIMULATION PARAMETERS ****************************/
#define XDIM 100
#define YDIM 100
#define DX 2.2649e9/XDIM
#define DY 2.2649e9/YDIM
#define MAX_ITERATIONS 5 //number of iterations to simulate (negative defers to MAX_TIME)
#define MAX_TIME -1.0 //upper bound on simulation time (negative defers to MAX_ITERATIONS)
/****************************************************************************/

/***************************** PHYSICS SETTINGS ****************************/
#define RADIATIVE_LOSSES_ON 1
#define AMBIENT_HEATING_ON 1
#define THERMAL_CONDUCTION_ON 1
#define FLUX_SATURATION 1 //caps thermal conductive flux at physical limit of plasma
/****************************************************************************/

/**************************** PHYSICAL PARAMETERS ***************************/
#define TEMP_CHROMOSPHERE 3.0e4 //chromospheric temperature, K (for radiation purposes)
#define RADIATION_RAMP 1.0e3 //width of falloff for low-temperature radiation, K
#define B0 100.0 //strength of B field at base of domain, Gauss
#define HEATING_RATE 1.0e-4 //uniform volumetric heating rate, erg cm^-3 s^-1
/****************************************************************************/

/*********************** TIME EVOLUTION SAFETY FACTORS **********************/
#define EPSILON 0.1 //time stepping safety factor (for CFL condition)
#define EPSILON_THERMAL 0.1 //safety factor for thermal conduction timestep (<0.5)
#define EPSILON_VISCOUS 1.0 //controls strength of artificial viscosity
#define EPSILON_RADIATIVE 0.1 //max fraction of total energy allowed to be lost in single radiative cycle
#define DT_THERMAL_MIN 1.0e-4 //minimum timestep for thermal conduction, s
#define RHO_MIN 1.0e-30 //min value for density (to avoid negative values)
#define TEMP_MIN 1.0e4 //minimum allowed temperature
#define THERMALENERGYFLOOR 1.0e-6 //min (nonzero) value for thermal energy density, erg cm^-3
/****************************************************************************/

/*************************** OUTPUT FILE SETTINGS ***************************/
#define OUTPUT_INTERVAL -1 //time steps between file outputs (set as negative to defer to TIME_OUTPUT_INTERVAL)
#define TIME_OUTPUT_INTERVAL 1.0 //simulation time between file outputs (set as negative to defer to OUTPUT_INTERVAL)
//Control which variables are output in .out file
#define RHO_OUT 1
#define TEMP_OUT 1
#define PRESS_OUT 1
#define RAD_OUT 1
#define ENERGY_OUT 1
#define VX_OUT 1
#define VY_OUT 1
#define DT_OUT 0
#define DT_THERMAL_OUT 1
#define DT_RAD_OUT 0
/****************************************************************************/

/*************************** BOUNDARY CONDITIONS ***************************/
//Boundary condition labels
#define PERIODIC 0
#define WALL 1
#define OPEN 2
//Boundary condition settings
#define XBOUND1 0
#define XBOUND2 0
#define YBOUND1 1
#define YBOUND2 2
/****************************************************************************/

/**************************** PHYSICAL CONSTANTS ****************************/
#define K_B 1.3807e-16 //boltzmann constant, erg K^-1
#define M_I 1.6726e-24 //ion mass, g
#define M_ELECTRON 9.1094e-28 //electron mass, g
#define BASE_GRAV 2.748e4 //acceleration due to gravity at surface, cm sec^-2
#define R_SUN 6.957e10 //radius of sun, cm
#define M_SUN 1.989e33 //mass of sun, g
#define GRAV_CONST 6.674e-8 //gravitational constant, dyn cm^2 g^-2
#define GAMMA 1.666667 //adiabatic index
#define KAPPA_0 1.0e-6 //thermal conductivity coefficient
#define PI 3.14159265358979323846
/****************************************************************************/

//Performance Benchmarking (Handle with care!)
//When turned on, outputs a .json file readable by chrome://tracing
//on Google Chrome for visual profiling purposes
#define BENCHMARKING_ON 0

#endif
