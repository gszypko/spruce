#include "antihelmholtz.hpp"

AntiHelmholtz::AntiHelmholtz(long double a, long double sep, long double dBdr,CurrentLoop::Axes ax):
    m_a{a}, m_sep{sep}, m_dBdr{dBdr}, m_axsym{ax}, m_conv{1.}
{
    // initialize loops
    m_loopR = CurrentLoop(m_a,-1,m_axsym,m_sep/2.);
    m_loopL = CurrentLoop(m_a,1,m_axsym,-m_sep/2.);

    // set linear magnetic field gradient along symmetry axis near the center of the anti-Helmholtz system
    std::vector<long double> r(3);
    r[m_axsym] = 1;
    std::vector<long double> B = get_field(r);
    m_conv = abs(m_dBdr/B[m_axsym]);
}

// return quadrupole magnetic field at given positioni
std::vector<long double> AntiHelmholtz::get_field(std::vector<long double> pos) const
{
    if (pos.size() != 3) std::cerr << "Position vector of format [x y z]" << std::endl;
    std::vector<long double> B1(3), B2(3), B(3);
    B1 = m_loopR.get_field(pos);
    B2 = m_loopL.get_field(pos);
    for (int i = 0; i < B.size(); i++){
        B[i] = (B1[i] + B2[i])*m_conv;
    }
    return B;
}