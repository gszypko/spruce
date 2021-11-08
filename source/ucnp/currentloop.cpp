#include "currentloop.hpp"

// constructor
CurrentLoop::CurrentLoop(long double a, long double I, Axes ax, long double pos):
    m_a{a}, m_I{I}, m_axsym{ax}, m_C{I/PI}, m_pos{pos}
{
    if (m_axsym == ax_z) m_axind = {0,1,2};
    else if (m_axsym == ax_y) m_axind = {2,0,1};
    else if (m_axsym == ax_x) m_axind = {1,2,0};
}

// copy constructor
CurrentLoop& CurrentLoop::operator=(const CurrentLoop& other)
{
    m_a = other.m_a;
    m_I = other.m_I;
    m_axsym = other.m_axsym;
    m_C = other.m_C;
    m_pos = other.m_pos;
    m_axind = other.m_axind;
    return *this;
}

// return magnetic field (in G) at position x, y, z (in cm)
std::vector<long double> CurrentLoop::get_field(std::vector<long double> pos) const
{
    // ensure length of position vector is 3
    if (pos.size() != 3) std::cerr << "Position vector of format [x y z]" << std::endl;

    // process input position depending on symmetry axis
    // the equations in this function assume z is the symmetry axis.
    // axind reorients the input position axes appropriately depending on the symmetry axis.
    // see the constructor for how the association is done
    pos[m_axsym] -= m_pos;
    long double x = pos[m_axind[0]];
    long double y = pos[m_axind[1]];
    long double z = pos[m_axind[2]];

    // define variables for getting magnetic field
    long double rho = sqrt(x*x+y*y);
    long double r = sqrt(x*x+y*y+z*z);
    long double alpha = sqrt(m_a*m_a+r*r-2.*m_a*rho);
    long double beta = sqrt(m_a*m_a+r*r+2.*m_a*rho);
    long double k = sqrt(1.-pow(alpha,2.)/pow(beta,2.));
    long double gamma = x*x-y*y;

    // get elliptic integrals
    long double eK = std::comp_ellint_1(k);
    long double eE = std::comp_ellint_2(k);

    // compute magnetic field
    std::vector<long double> B(3);
    long double fac = m_C/(2.*alpha*alpha*beta);
    B[m_axind[0]] = fac*x*z/(rho*rho)*((m_a*m_a+r*r)*eE-alpha*alpha*eK);
    B[m_axind[1]] = fac*y*z/(rho*rho)*((m_a*m_a+r*r)*eE-alpha*alpha*eK);
    B[m_axind[2]] = fac*((m_a*m_a-r*r)*eE+alpha*alpha*eK);

    // handle singularity on the symmetry axis by using equations for near-axis approximation
    if (rho/m_a < 1e-5){
        B[m_axind[0]] = 3*PI*m_a*m_a*x*z/(4*(m_a*m_a+z*z));
        B[m_axind[1]] = 3*PI*m_a*m_a*y*z/(4*(m_a*m_a+z*z));
    }

    return B;
}