#include "ucnputils.hpp"
#include "constants.hpp"
#include <cmath>

namespace phys
{
// ratio between average Coulomb interaction energy and thermal energy (inputs in cgs)
double coulomb_coupling(double n,double T)
// n: plasma density (cm^-3)
// T: temperature (K) of a particular species (i.e., electron or ion)
{ return nearest_coulomb_pot(n)/(K_B*T); }

// average interparticle spacing of gas with density n (SI units)
double wigner_seitz_radius(double n) 
// n: plasma density (cm^-3)
{ return pow(3./(4.*PI*n),1./3.); } 

// plasma oscillation frequency (rad/s, inputs cgs)
double plasma_freq(double n,double m) 
// n: plasma density (cm^-3)
// m: particle mass (g)
{ return sqrt(4*PI*n*pow(E,2.)/m); }

// einsten frequency (rad/s, inputs cgs)
double einstein_freq(double n,double m)
// n: plasma density (cm^-3)
// m: particle mass (g)
{ return plasma_freq(n,m)/sqrt(3.); }

// average Coulomb potential between nearest neighbors, cgs units
double nearest_coulomb_pot(double n) 
// n: plasma density (cm^-3)
{ return E*E/wigner_seitz_radius(n); } 

// debye length of a particular plasma species, cgs units
double debye_length(double n,double T) 
// n: plasma density (cm^-3)
// T: temperature (K) of a particular plasma species
{ return sqrt(K_B*T/(4*PI*n*E*E));}

// plasma screening parameter for electrons, dimensionless, inputs cgs
double screening_parameter(double n, double Te)
// n: plasma density (cm^-3)
// Te: electron temperature (K)
{ return wigner_seitz_radius(n)/debye_length(n,Te); } 

// equilibrium plasma temperature due to disorder-induced heating, inputs cgs
double dih_temp(double n,double Te) 
// n: plasma density (cm^-3)
// Te: electron temperature (K)
{ return 2.*nearest_coulomb_pot(n)/3./K_B*(1.+screening_parameter(n,Te)/2.); }

// timescale for a Gaussian UCNP expansion into vacuum, cgs units
double tau_exp(double sig,double m, double Te)
// n: plasma density (cm^-3)
// m: ion mass (g)
// Te: electron temperature (K)
{
    return sqrt(m*pow(sig,2.)/(K_B*Te));
}
}