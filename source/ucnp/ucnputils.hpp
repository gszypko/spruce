#ifndef UCNPUTILS_HPP
#define UCNPUTILS_HPP

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