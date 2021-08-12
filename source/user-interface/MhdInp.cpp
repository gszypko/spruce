#include "MhdInp.hpp"

// class constructor
MhdInp::MhdInp(size_t Nx,size_t Ny)  
{
    m_Nx = Nx; m_Ny = Ny;
    m_grids.resize(num_inp,Grid(m_Nx,m_Ny,0.));
    m_initialized.resize(num_inp,false);
}

// set an element of <m_grids> corresponding to integer input
void MhdInp::set_var(const int& var,const Grid& grid)
{
    if (grid.rows() != m_Nx || grid.cols() != m_Ny) std::cerr << "Input grid dimensions for m_grid[" << var << "]." << std::endl;
    m_grids[var] = grid;
    m_initialized[var] = true;
}

// check if all elements of <m_grids> were initialized
void MhdInp::all_initialized() const
{
    for (int i = 0; i < num_inp; i++){
        if (!m_initialized[i]) std::cerr << "The grid for <" << m_varnames[i] << "> was not initialized." << std::endl;
    } 
}

// return all of the input grids, but first check that they were all initialized
std::vector<Grid> MhdInp::grids()
{
    all_initialized();
    return m_grids;
}