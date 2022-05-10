#include "MhdInp.hpp"
#include "plasmadomain.hpp"
#include <stdexcept>
#include <cctype>

// class constructor
MhdInp::MhdInp(size_t Nx,size_t Ny)  
{
    m_Nx = Nx; m_Ny = Ny;
    for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_grids.push_back(Grid(m_Nx,m_Ny,0.0));
    for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_initialized.push_back(false);
    m_ion_mass = 0.0; m_adiabatic_index = 0.0; m_duration = 0.0; m_time = 0.0;
}

// Copy assignment operator
MhdInp& MhdInp::operator=(const MhdInp& other)
{
    m_Nx = other.m_Nx; m_Ny = other.m_Ny;
    m_adiabatic_index = other.m_adiabatic_index; m_ion_mass = other.m_ion_mass; 
    m_duration = other.m_duration; m_time = other.m_time;
    for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_grids.push_back(Grid(m_Nx,m_Ny,0.0));
    for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_initialized.push_back(false);

    for(int i : PlasmaDomain::state_vars){
        m_grids[i] = other.m_grids[i];
        m_initialized[i] = other.m_initialized[i];
    }

    return *this;
}

// set an element of <m_grids> corresponding to integer input
// When setting d_x or d_y, origin_pos must be specified as "lower", "center", or "upper"
// to define to location of the origin in the domain for the generated pos_x and pos_y
void MhdInp::set_var(int var,const Grid& grid,const std::string origin_pos)
{
    if (grid.rows() != m_Nx || grid.cols() != m_Ny) std::cerr << "Input grid dimensions for m_grid[" << var << "]." << std::endl;
    m_grids[var] = grid;
    m_initialized[var] = true;
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

void MhdInp::set_time(double time)
{
    m_time = time;
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

double MhdInp::time()
{
    return m_time;
}

// write grids to state file
void MhdInp::write_state_file(const fs::path& directory) const
{
    all_initialized();
    if (!fs::exists(directory)) fs::create_directories(directory);
    fs::path filename = "init.state";
    std::ofstream outfile(directory/filename);
    std::string dlm = ",";
    outfile << "xdim,ydim" << std::endl;
    outfile << m_Nx << dlm << m_Ny << std::endl;
    outfile << "ion_mass" << std::endl;
    outfile << m_ion_mass << std::endl;
    outfile << "adiabatic_index" << std::endl;
    outfile << m_adiabatic_index << std::endl;
    outfile << "t=" << m_time << std::endl;
    outfile << "duration=" << m_duration << std::endl;
    for (auto& i : PlasmaDomain::state_vars){
        outfile << PlasmaDomain::m_var_names[i] << std::endl;
        outfile << m_grids[i].format(',','\n',-1);
    }
}

// check if all elements of <m_grids> were initialized
void MhdInp::all_initialized() const
{
    for (int i : PlasmaDomain::state_vars){
        if (!m_initialized[i]){
            std::cerr << "The grid for <" << PlasmaDomain::m_var_names[i] << "> was not initialized." << std::endl;
            exit(1);
        }
    } 
}

// return all of the input grids, but first check that they were all initialized
std::vector<Grid> MhdInp::grids()
{
    all_initialized();
    return m_grids;
}