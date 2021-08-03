#ifndef CONSTANTS_UCNP_H
#define CONSTANTS_UCNP_H

// constants in SI units
namespace cts
{
    // general constants
    constexpr double kB {1.3806e-23};       // Boltzmann constant (units: J/K)
    constexpr double eps0 {8.8542e-12};     // vacuum permittivity (units: F/m)
    constexpr double u0 {1.2566e-6};        // vacuum permeability (units: H/m)
    constexpr double uB {9.274e-24};        // Bohr magneton (units: J/T)
    constexpr double h {6.6261e-34};        // Planck's constant (units: Js)
    constexpr double hbar {1.0546e-34};     // hbar (units: Js)
    constexpr double c {299792458};         // speed of light (units: m/s)
    constexpr double pi {3.14159265359};    // pi

    // electron properties
    constexpr double mE {9.1095e-31};   // electron mass (units: kg)
    constexpr double e {1.602e-19};     // electron charge (units: C)

    // Sr+ properties
    constexpr double mI {1.455e-25};    // Sr+ mass (units: kg)

    // Sr+ transition natural linewidths (units: s^-1)
    constexpr double gam422 {1.2755e8};     // natural linewidth of 2S_1/2 -> 2P_1/2 422 nm imaging transition
    constexpr double gam408 {1.4074e8};     // natural linewdith of 2S_1/2 -> 2P_3/2 408 nm cooling transition
    constexpr double gam1033 {8.6959e6};    // natural linewidth of 2D_5/2 -> 2P_3/2 1033 nm repump transition
    constexpr double gam1092 {7.46e6};      // natural linewidth of 2D_3/2 -> 2P_1/2 1092 nm repump transition

    // 422 nm laser linewidth (units: s^-1)
    constexpr double gamL {3.1416e7};

    // Sr+ transition wavelengths (units: m)
    constexpr double lam422 {421.552e-9};      // wavelength of imaging transition
    constexpr double lam408 {407.771e-9};      // wavelength of cooling transition
    constexpr double lam1033 {1032.7311e-9};   // wavelength of strong repump transition
    constexpr double lam1092 {1091.4887e-9};   // wavelength of weak repump transition
}

#endif
