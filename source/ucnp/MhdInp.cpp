#include "MhdInp.hpp"

// class constructor
MhdInp::MhdInp(size_t Nx,size_t Ny,PlasmaDomain& pd,std::string eqs_set):
    m_Nx{Nx}, m_Ny{Ny}, m_eqs{EquationSet::spawnEquationSet(pd,eqs_set)}
{  
    // collect full list of variable names - including plasma domain and equation set grids  
    for (const auto& name : PlasmaDomain::m_internal_var_names) m_grid_names.push_back(name);
    for (const auto& ind : m_eqs->state_variables()) m_grid_names.push_back(m_eqs->nameFromIndex(ind));    
    // initialize grids for each variable and indicate they have not yet been initialized
    for (int i=0; i<m_grid_names.size(); i++){
        m_grid_indices[m_grid_names[i]] = i;
        m_grids.push_back(Grid::Zero(1,1));
        m_initialized.push_back(false);
    }
}

int MhdInp::name2index(std::string name) const
{
    int ind{};
    try { 
        ind = m_grid_indices.at(name); }
    catch (const std::out_of_range& e) {
        std::cerr << "Variable name " << name << " not recognized" << std::endl;
        assert(false);
    }
    return ind;
}

// Copy assignment operator
// MhdInp& MhdInp::operator=(const MhdInp& other)
// {
//     m_Nx = other.m_Nx; m_Ny = other.m_Ny;
//     m_adiabatic_index = other.m_adiabatic_index; m_ion_mass = other.m_ion_mass; 
//     m_duration = other.m_duration; m_time = other.m_time;
//     for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_grids.push_back(Grid(m_Nx,m_Ny,0.0));
//     for(int i=0; i<PlasmaDomain::state_vars.size(); i++) m_initialized.push_back(false);

//     for(int i : PlasmaDomain::state_vars){
//         m_grids[i] = other.m_grids[i];
//         m_initialized[i] = other.m_initialized[i];
//     }

//     return *this;
// }

// set an element of <m_grids> corresponding to integer input
// When setting d_x or d_y, origin_pos must be specified as "lower", "center", or "upper"
// to define to location of the origin in the domain for the generated pos_x and pos_y
void MhdInp::set_var(std::string grid_name,const Grid& grid)
{
    // ensure that grid dimensions are compatible
    if (grid.rows() != m_Nx || grid.cols() != m_Ny){
        std::cerr << "Input grid dimensions for m_grid[" << grid_name << "]." << std::endl;
        assert(false);
    } 

    int ind {name2index(grid_name)};
    m_grids[ind] = grid;
    m_initialized[ind] = true;
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

void MhdInp::set_time_output_interval(double time_output_interval)
{
    m_time_output_interval = time_output_interval;
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
    for (int i=0; i<m_grids.size(); i++){
        outfile << m_grid_names[i] << std::endl;
        outfile << m_grids[i].format(',','\n',-1);
    }
}

void MhdInp::read_config(const fs::path& filepath)
{
    assert(filepath.extension().string()==".config");
    assert(fs::exists(filepath) && !fs::is_empty(filepath) && "File must exist and not be empty.");
    m_config_data.reserve(100);
    std::ifstream fileStream(filepath);
    while (fileStream.good()){
        std::string currLine;
        std::getline(fileStream,currLine);
        m_config_data.push_back(currLine);
        if (m_config_data.size() == m_config_data.capacity()) 
            m_config_data.reserve(2*m_config_data.capacity());
    }
    fileStream.close();
    m_config_data.shrink_to_fit();
}

void MhdInp::update_config(std::string name,std::string val)
{
    // ensure that input is a config name
    auto it = std::find(m_config_names.begin(),m_config_names.end(),name);
    assert(it!=m_config_names.end());
    // find the line in config data that starts with this value
    bool config_found{false};
    for (int i=0; i<m_config_data.size(); i++){
        std::string curr_line = m_config_data[i];
        std::string lhs,rhs;
        std::istringstream ss(curr_line);
        std::getline(ss,lhs,'=');
        std::getline(ss,rhs);
        std::string lhs_cleared = lhs; lhs_cleared.erase(std::remove(lhs_cleared.begin(),lhs_cleared.end(),' '),lhs_cleared.end());
        if (lhs_cleared.compare(name)==0){
             m_config_data[i] = lhs + " = " + val;
             config_found = true;
        }
    }
    if (!config_found) m_config_data.push_back(name + " = " + val);
    
}

void MhdInp::write_config_file(const fs::path& directory) const
{
    if (!fs::exists(directory)) fs::create_directories(directory);
    fs::path filename = "ucnp.config";
    std::ofstream outfile(directory/filename);
    for (const auto& curr_line : m_config_data){
        outfile << curr_line << std::endl;
    }
    outfile.close();
}

// check if all elements of <m_grids> were initialized
void MhdInp::all_initialized() const
{
    for (int i=0; i<m_grids.size(); i++){
        if (!m_initialized[i]){
            std::cerr << "The grid for <" << m_grid_names[i] << "> was not initialized." << std::endl;
            assert(false);
        }
    } 
}

bool MhdInp::is_config(const std::string& name) const
{
    auto it = std::find(m_config_names.begin(),m_config_names.end(),name);
    return it != m_config_names.end();
}

// return all of the input grids, but first check that they were all initialized
std::vector<Grid> MhdInp::grids()
{
    all_initialized();
    return m_grids;
}