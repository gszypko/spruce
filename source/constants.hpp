#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#define N_GHOST 2 //number of ghost cells to use for Open and Wall boundaries

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
#define M_SR 1.455e-28 // mass of sr ion, g
#define C 29979245800 // speed of light, cm/s
#define E 4.80320425e-10 // electric charge, CGS: statcoulombs
/****************************************************************************/

//Performance Benchmarking (Handle with care!)
//When turned on, outputs a .json file readable by chrome://tracing
//on Google Chrome for visual profiling purposes
#define BENCHMARKING_ON 0

#endif
