#include "MhdInp.hpp"
#include "plasmadomain.hpp"

// class constructor
MhdInp::MhdInp(size_t Nx,size_t Ny)  
{
    m_Nx = Nx; m_Ny = Ny;
    m_grids.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,Grid(m_Nx,m_Ny,0.));
    m_initialized.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,false);
    m_ion_mass = 0.0; m_adiabatic_index = 0.0; m_duration = 0.0;
}

// Copy assignment operator
MhdInp& MhdInp::operator=(const MhdInp& other)
{
    m_Nx = other.m_Nx; m_Ny = other.m_Ny;
    m_adiabatic_index = other.m_adiabatic_index; m_ion_mass = other.m_ion_mass; m_duration = other.m_duration;
    m_grids.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,Grid(m_Nx,m_Ny,0.));
    m_initialized.resize(PlasmaDomain::state_var_end - PlasmaDomain::state_var_start,false);

    for(int i=PlasmaDomain::state_var_start; i<PlasmaDomain::state_var_end; i++){
        m_grids[i] = other.m_grids[i];
        m_initialized[i] = other.m_initialized[i];
    }

    return *this;
}

// set an element of <m_grids> corresponding to integer input
// When setting d_x or d_y, origin_pos must be specified as "lower", "center", or "upper"
// to define to location of the origin in the domain for the generated pos_x and pos_y
void MhdInp::set_var(const int& var,const Grid& grid,const std::string origin_pos)
{
    if (grid.rows() != m_Nx || grid.cols() != m_Ny) std::cerr << "Input grid dimensions for m_grid[" << var << "]." << std::endl;
    if(var == PlasmaDomain::pos_x || var == PlasmaDomain::pos_y){
        std::cerr << "Position grids are only to be set by way of defining d_x and d_y" << std::endl;
        return;
    }
    m_grids[var] = grid;
    m_initialized[var] = true;
    if(var == PlasmaDomain::d_x){
        if(origin_pos == "") std::cerr << "Warning: x origin position not specified; default is lower bound" << std::endl;
        m_grids[PlasmaDomain::pos_x] = PlasmaDomain::convertCellSizesToCellPositions(grid,0,origin_pos);
        m_initialized[PlasmaDomain::pos_x] = true;
    }
    else if(var == PlasmaDomain::d_y){
        if(origin_pos == "") std::cerr << "Warning: y origin position not specified; default is lower bound" << std::endl;
        m_grids[PlasmaDomain::pos_y] = PlasmaDomain::convertCellSizesToCellPositions(grid,1,origin_pos);
        m_initialized[PlasmaDomain::pos_y] = true;
    }
}

void MhdInp::set_ion_mass(double mass)
{
    m_ion_mass = mass;
}

void MhdInp::set_adiabatic_index(double index)
{
    m_adiabatic_index = index;
}

void MhdInp::set_duration(double duration)
{
    m_duration = duration;
}

double MhdInp::ion_mass()
{
    if(m_ion_mass == 0.0) std::cerr << "Ion mass was not initialized." << std::endl;
    return m_ion_mass;
}

double MhdInp::adiabatic_index()
{
    if(m_adiabatic_index == 0.0) std::cerr << "Adiabatic index was not initialized." << std::endl;
    return m_adiabatic_index;
}

double MhdInp::duration()
{
    if(m_duration == 0.0) std::cerr << "Duration was not initialized." << std::endl;
    return m_duration;
}

// check if all elements of <m_grids> were initialized
void MhdInp::all_initialized() const
{
    for (int i = PlasmaDomain::state_var_start; i < PlasmaDomain::state_var_end; i++){
        if (!m_initialized[i]) std::cerr << "The grid at index " << i << " was not initialized." << std::endl;
    } 
}

// return all of the input grids, but first check that they were all initialized
std::vector<Grid> MhdInp::grids()
{
    all_initialized();
    return m_grids;
}