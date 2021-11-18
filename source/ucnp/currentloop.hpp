#ifndef CURRENTLOOP_HPP
#define CURRENTLOOP_HPP

#include "constants.hpp"
#include <iostream> // std::cerr, std::endl, std::cout
#include <math.h> // std::sqrt, std::abs
#include <vector> // std::vector
#include <string> // std::string, std::getline
#include <filesystem> // std::filesystem
namespace fs = std::filesystem;
#include <fstream> // std::ofstream, std::ifstream
#include <sstream> // std::istringstream

// this class provides the ability to calculate the magnetic field produced by a circular current loop
// of given radius, symmetry axis, and position relative to the origin of the symmetry axis.

// all units are cgs

// equations taken from "Simple analytic expressions for the magnetic field of a circular current loop" by J. Simpson

class CurrentLoop
{
public: 
    enum Axes{ax_x,ax_y,ax_z,ax_num};   // axis information
    
    CurrentLoop() : CurrentLoop(1.,1.,Axes::ax_x,1.) {}
    CurrentLoop(long double a, long double I, Axes ax, long double pos);
    CurrentLoop& operator=(const CurrentLoop& other);

    std::vector<long double> get_field(std::vector<long double> pos) const;

private:
    // member variables to be initialized with constructor

    long double m_a; // radius of current loop in cm
    long double m_I; // current through loop in statampere
    long double m_C; // constant of proportionality for magnetic field
    Axes m_axsym; // symmetry axis for this loop
    long double m_pos; // position of coil along symmetry axis
    std::vector<long double> m_axind; // associates input position axes
};

#endif