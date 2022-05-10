#ifndef MHD_INP_HPP
#define MHD_INP_HPP

#include "grid.hpp"

#include <iostream> // std::cerr, std::endl, std::cout
#include <vector> // std::vector
#include <string> // std::string, std::getline
#include <filesystem>
#include <fstream>

namespace fs = std::filesystem;

class MhdInp
{
public:
    MhdInp& operator=(const MhdInp& other);
    MhdInp(size_t Nx,size_t Ny);
    MhdInp() : MhdInp(1,1) {}
    void set_var(int var,const Grid& grid,const std::string origin_pos = "");
    void set_ion_mass(double mass);
    void set_adiabatic_index(double index);
    void set_duration(double duration);
    void set_time(double time);
    void all_initialized() const;
    std::vector<Grid> grids();
    double ion_mass();
    double adiabatic_index();
    double duration();
    double time();
    void write_state_file(const fs::path& directory) const;
private:
    size_t m_Nx; // number of elements along x-axis contained in <m_grids> (i.e., the rows of the grid)
    size_t m_Ny; // number of elements along y-axis contained in <m_grids> (i.e., the columns of the grid)
    std::vector<Grid> m_grids; // input grids for MHD simulation - enumed by <Vars>
    std::vector<bool> m_initialized; // records whether the user has initialized each element of <m_grids>
    double m_ion_mass;
    double m_adiabatic_index;
    double m_duration;
    double m_time;
};

#endif