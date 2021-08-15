#include "MhdInp.hpp"
#include "plasmadomain.hpp"

// class constructor
MhdInp::MhdInp(size_t Nx,size_t Ny)  
{
    m_Nx = Nx; m_Ny = Ny;
    m_grids.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,Grid(m_Nx,m_Ny,0.));
    m_initialized.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,false);
}

// Copy assignment operator
MhdInp& MhdInp::operator=(const MhdInp& other)
{
    m_Nx = other.m_Nx; m_Ny = other.m_Ny;
    m_grids.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,Grid(m_Nx,m_Ny,0.));
    m_initialized.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,false);

    for(int i=PlasmaDomain::state_var_start; i<PlasmaDomain::state_var_end; i++){
        m_grids[i] = other.m_grids[i];
        m_initialized[i] = other.m_initialized[i];
    }

    return *this;
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
    for (int i = PlasmaDomain::state_var_start; i < PlasmaDomain::state_var_end; i++){
        if (!m_initialized[i]) std::cerr << "The grid for <" << m_varnames[i] << "> was not initialized." << std::endl;
    } 
}

// return all of the input grids, but first check that they were all initialized
std::vector<Grid> MhdInp::grids()
{
    all_initialized();
    return m_grids;
}