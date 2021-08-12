#ifndef MHD_INP_HPP
#define MHD_INP_HPP

#include "grid.hpp"

#include <iostream> // std::cerr, std::endl, std::cout
#include <vector> // std::vector
#include <string> // std::string, std::getline

class MhdInp
{
public:
    enum Vars {x,y,rho,temp,mom_x,mom_y,b_x,b_y,b_z,num_inp}; // enum for <m_grids> and <m_initialized>
    MhdInp& operator=(const MhdInp& other);
    MhdInp(size_t Nx,size_t Ny);
    void set_var(const int& var,const Grid& grid);
    void all_initialized() const;
    std::vector<Grid> grids();
private:
    size_t m_Nx; // number of elements along x-axis contained in <m_grids>
    size_t m_Ny; // number of elements along y-axis contained in <m_grids>
    const std::vector<std::string> m_varnames { "x","y","rho","temp","mom_x","mom_y","b_x","b_y","b_z" };
    std::vector<Grid> m_grids; // input grids for MHD simulation - enumed by <Vars>
    std::vector<bool> m_initialized; // records whether the user has initialized each element of <m_grids>
};

#endif