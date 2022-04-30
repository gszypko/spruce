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
    m_ion_mass = 0.0; m_adiabatic_index = 0.0; m_duration = 0.0;
}

// Copy assignment operator
MhdInp& MhdInp::operator=(const MhdInp& other)
{
    m_Nx = other.m_Nx; m_Ny = other.m_Ny;
    m_adiabatic_index = other.m_adiabatic_index; m_ion_mass = other.m_ion_mass; m_duration = other.m_duration;
    for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_grids.push_back(Grid(m_Nx,m_Ny,0.0));
    for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_initialized.push_back(false);

    for(int i : PlasmaDomain::state_vars){
        m_grids[i] = other.m_grids[i];
        m_initialized[i] = other.m_initialized[i];
    }

    return *this;
}

// void MhdInp::read_custom_input_csvs(const fs::path &grid_path)
// {
//     for(auto const& dir_entry : fs::directory_iterator{grid_path}){
//         fs::path file_path = dir_entry.path();
//         std::string exten = file_path.extension().string();
//         if(exten == ".csv"){
//             std::string varname = file_path.stem().string();
//             //Ensure that variable in file is valid
//             auto it = std::find(PlasmaDomain::m_var_names.begin(),PlasmaDomain::m_var_names.end(),varname);
//             if(it != PlasmaDomain::m_var_names.end()){
//                 auto index = std::distance(PlasmaDomain::m_var_names.begin(),it);
//                 //Read in Grid corresponding to variable
//                 std::ifstream in_file(file_path.string());
//                 std::string row; std::string el;
//                 std::vector<double> data_vec(0);
//                 int i=0;
//                 while(std::getline(in_file,row)){
//                     std::istringstream ss_row(row);
//                     int j=0;
//                     double curr_val;
//                     while(std::getline(ss_row,el,',')){
//                         try { curr_val = std::stod(el); }
//                         catch (const std::invalid_argument&){
//                             std::cerr << "Issue reading " << file_path.filename().string() << ", ";
//                             std::cerr << "possibly due to BOMs or other encoding characters.\n";
//                             std::cerr << "Attempting to repair... ";
//                             bool success = false;
//                             while (!success && el.size() > 1){
//                                 success = true;
//                                 el.erase(0,1); 
//                                 try { curr_val = std::stod(el); } 
//                                 catch(const std::invalid_argument&) { 
//                                     success = false;
//                                 }
//                             }
//                             if (!success){
//                                 std::cerr << "Repair unsuccessful. Try removing any encoding characters from " << file_path.filename().string() << std::endl;
//                                 throw;
//                             }
//                             else std::cerr << "Repair successful.\n";
//                         }
//                         data_vec.push_back(curr_val);
//                         j++;
//                     }
//                     if(m_Ny == 1) m_Ny = j;
//                     else assert(j == m_Ny && "All rows must have the same length in all custom input files (i.e. consistent y-dimension)");
//                     i++;
//                 }
//                 if(m_Nx == 1) m_Nx = i;
//                 else assert(i == m_Nx && "All custom input files must have the same number of rows (i.e. consistent x-dimension)");
//                 Grid curr_grid(m_Nx,m_Ny,data_vec);
//                 set_var(index,curr_grid);
//             } else std::cout << "Variable name " << varname << " not recognized in custom input parsing\n";
//         }
//     }
// }

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
    for (int i : PlasmaDomain::state_vars){
        if (!m_initialized[i]) std::cerr << "The grid at index " << i << " was not initialized." << std::endl;
    } 
}

// return all of the input grids, but first check that they were all initialized
std::vector<Grid> MhdInp::grids()
{
    all_initialized();
    return m_grids;
}