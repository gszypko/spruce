#ifndef ANTIHELMHOLTZ_HPP
#define ANTIHELMHOLTZ_HPP

#include "currentloop.hpp"

// compute magnetic field for two-coil anti-helmholtz configuration. all units are cgs
class AntiHelmholtz
{
public:    
    AntiHelmholtz(long double a, long double sep, long double dBdr, CurrentLoop::Axes ax);
    
    std::vector<long double> get_field(std::vector<long double> pos) const;

private:
    // member variables to be initialized with constructor

    long double m_a; // radius of current loop in cm
    long double m_sep; // separation between two coils
    long double m_dBdr; // linear magnetic field gradient along symmetry axis (near center)
    CurrentLoop::Axes m_axsym; // symmetry axis for this loop
    long double m_conv; 

    CurrentLoop m_loopR;
    CurrentLoop m_loopL;
};

#endif