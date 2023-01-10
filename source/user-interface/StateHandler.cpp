#include "StateHandler.hpp"
#include <cassert>
#include <fstream>
#include <cmath>
#include "antihelmholtz.hpp"

StateHandler::StateHandler(const std::string& eqn_set_name):
    m_eqn_set_name{eqn_set_name}
{
    m_eqs = EquationSet::instantiateDefault(m_pd,m_eqn_set_name);
    for (const auto& name : PlasmaDomain::m_gridnames) m_gridnames.push_back(name);
    for (const auto& ind : m_eqs->state_variables()) m_gridnames.push_back(m_eqs->index2name(ind));
    m_grids.resize(m_gridnames.size(),Grid::Zero(1,1));
    m_grids_initialized.resize(m_gridnames.size(),false);
}

int StateHandler::gridname2ind(const std::string& name) const
{
    auto it = std::find(m_gridnames.begin(),m_gridnames.end(),name);
    if (it == m_gridnames.end()){
        std::cerr << "<" << name << "> is not a valid grid name." << std::endl;
        assert(false);
    }
    return std::distance(m_gridnames.begin(),it);
}

int StateHandler::varname2ind(const std::string& name) const
{
    auto it = std::find(m_varnames.begin(),m_varnames.end(),name);
    if (it == m_varnames.end()){
        std::cerr << "<" << name << "> is not a valid var name." << std::endl;
        assert(false);
    }
    return std::distance(m_varnames.begin(),it);
}

void StateHandler::setvar(const std::string& name, double val)
{
    int ind = varname2ind(name);
    m_vars[ind] = val;
    m_vars_initialized[ind] = true;
}

void StateHandler::setgrid(const std::string& name, Grid grid)
{
    int ind = gridname2ind(name);
    m_grids[ind] = grid;
    m_grids_initialized[ind] = true;
}

double StateHandler::getvar(const std::string& name) const
{
    int ind = varname2ind(name);
    assert(m_vars_initialized[ind]);
    return m_vars[ind];
}

Grid StateHandler::getgrid(const std::string& name) const
{
    int ind = gridname2ind(name);
    assert(m_grids_initialized[ind]);
    return m_grids[ind];
}

void StateHandler::all_initialized() const
{
    for (int i=0; i<m_varnames.size(); i++){
        if (m_vars_initialized[i] == false){
            std::cerr << "Varible <" << m_varnames[i] << "> is not initialized." << std::endl;
            std::cerr << "All grids and vars must be initialized before writing state file." << std::endl;
            assert(false);
        }
    }
    for (int i=0; i<m_gridnames.size(); i++){
        if (m_grids_initialized[i] == false){
            std::cerr << "Grid <" << m_gridnames[i] << "> is not initialized." << std::endl;
            std::cerr << "All grids and vars must be initialized before writing state file." << std::endl;
            assert(false);
        }
    }
}

void StateHandler::write_state_file(const fs::path& directory) const
{
    all_initialized();
    if (!fs::exists(directory)) fs::create_directories(directory);
    fs::path filename = "init.state";
    std::ofstream outfile(directory/filename);
    outfile << "xdim,ydim" << std::endl;
    outfile << getvar("xdim") << "," << getvar("ydim") << std::endl;
    outfile << "ion_mass" << std::endl;
    outfile << getvar("ion_mass") << std::endl;
    outfile << "adiabatic_index" << std::endl;
    outfile << getvar("adiabatic_index") << std::endl;
    outfile << "t=" << getvar("time") << std::endl;
    for (int i=0; i<m_grids.size(); i++){
        outfile << m_gridnames[i] << std::endl;
        outfile << m_grids[i].format(',','\n',-1);
    }
}

void StateHandler::setup(const std::unique_ptr<Settings>& pms)
{
    //*** initialize variables
    setvar("xdim",pms->getval("Nx"));
    setvar("ydim",pms->getval("Ny"));
    setvar("ion_mass",pms->getval("m_i"));
    setvar("adiabatic_index",pms->getval("gam"));
    setvar("time",0.);
    //*** initialize grids for plasma domain
    // position grids
    Grid d_x, d_y, pos_x, pos_y;
    Grid dx_vec = Grid(1,getvar("xdim"),2*pms->getval("x_lim")/(getvar("xdim")-1));
    Grid dy_vec = Grid(1,getvar("ydim"),2*pms->getval("y_lim")/(getvar("ydim")-1));
    Grid::MeshGrid(dx_vec,dy_vec,d_x,d_y);
    pos_x = PlasmaDomain::convertCellSizesToCellPositions(d_x,0,"center");
    pos_y = PlasmaDomain::convertCellSizesToCellPositions(d_y,1,"center");   
    setgrid("d_x",d_x);
    setgrid("d_y",d_y);
    setgrid("pos_x",pos_x);
    setgrid("pos_y",pos_y);
    // external magnetic fields
    AntiHelmholtz quad(30.,120.,pms->getval("dBdx"),CurrentLoop::ax_x);
    std::vector<Grid> B(3.);
    for (int i = 0; i < B.size(); i++){ // for each B-field component (x, y, z)
        B[i] = Grid::Zero(getvar("xdim"),getvar("ydim")); // initialize that component with zeros
        for (int j = 0; j < B[i].rows(); j++){
            for (int k = 0; k < B[i].cols(); k++){
                std::vector<long double> Bquad = quad.get_field({pos_x(j,k),pos_y(j,k),0.});
                B[i](j,k) = Bquad[i];
            }
        }
    } 
    setgrid("be_x",B[0]);
    setgrid("be_y",B[1]);
    setgrid("be_z",B[2]);
    //*** initialize grids for equation set
    if (m_eqn_set_name == "ideal_mhd") setup_idealmhd(pms);
    else if (m_eqn_set_name == "ideal_mhd_cons") setup_idealmhdcons(pms);
    else if (m_eqn_set_name == "ideal_mhd_2E") setup_idealmhd2e(pms);
    else if (m_eqn_set_name == "ideal_2F") setup_ideal2F(pms);
    else assert(false && "Equation set not recognized.");
}

void StateHandler::setup_idealmhd(const std::unique_ptr<Settings>& pms)
{
    Grid n = setup_density(pms);
    Grid zeros = Grid::Zero(getvar("xdim"),getvar("ydim"));
    setgrid("rho",n*getvar("ion_mass"));
    setgrid("temp",Grid(getvar("xdim"),getvar("ydim"),pms->getval("Te")));
    setgrid("mom_x",zeros);
    setgrid("mom_y",zeros);
    setgrid("bi_x",zeros);
    setgrid("bi_y",zeros);
    setgrid("grav_x",zeros);
    setgrid("grav_y",zeros);
}

void StateHandler::setup_idealmhdcons(const std::unique_ptr<Settings>& pms)
{
    setup_idealmhd(pms);
}

void StateHandler::setup_idealmhd2e(const std::unique_ptr<Settings>& pms)
{
    Grid n = setup_density(pms);
    Grid zeros = Grid::Zero(getvar("xdim"),getvar("ydim"));
    setgrid("rho",n*getvar("ion_mass"));
    setgrid("e_temp",Grid(getvar("xdim"),getvar("ydim"),pms->getval("Te")));
    setgrid("i_temp",Grid(getvar("xdim"),getvar("ydim"),pms->getval("Ti")));
    setgrid("mom_x",zeros);
    setgrid("mom_y",zeros);
    setgrid("bi_x",zeros);
    setgrid("bi_y",zeros);
    setgrid("grav_x",zeros);
    setgrid("grav_y",zeros);
}

void StateHandler::setup_ideal2F(const std::unique_ptr<Settings>& pms)
{
    Grid n = setup_density(pms);
    Grid zeros = Grid::Zero(getvar("xdim"),getvar("ydim"));
    setgrid("i_rho",n*getvar("ion_mass"));
    setgrid("e_rho",n*M_ELECTRON);
    setgrid("i_mom_x",zeros);
    setgrid("i_mom_y",zeros);
    setgrid("e_mom_x",zeros);
    setgrid("e_mom_y",zeros);
    setgrid("i_temp",Grid(getvar("xdim"),getvar("ydim"),pms->getval("Ti")));
    setgrid("e_temp",Grid(getvar("xdim"),getvar("ydim"),pms->getval("Te")));
    setgrid("bi_x",zeros);
    setgrid("bi_y",zeros);
    setgrid("bi_z",zeros);
    setgrid("E_x",zeros);
    setgrid("E_y",zeros);
    setgrid("E_z",zeros);
    setgrid("grav_x",zeros);
    setgrid("grav_y",zeros);
}

Grid StateHandler::setup_density(const std::unique_ptr<Settings>& pms) const
{
    double n_max { pms->getval("n") };
    double n_min { pms->getval("n_min") };
    Grid n{Grid::Zero(getvar("xdim"),getvar("ydim"))};
    // choose overall density distribution for plasma
    if (pms->getopt("n_dist") == "gaussian") 
        n = Grid::Gaussian2D(getgrid("pos_x"),getgrid("pos_y"),n_max,n_min,pms->getval("sig_x"),pms->getval("sig_y"),0,0);
    else if (pms->getopt("n_dist")=="exponential") 
        n = Grid::Exp2D(getgrid("pos_x"),getgrid("pos_y"),n_max,n_min,pms->getval("sig_x"),pms->getval("sig_y"),0,0);
    else if (pms->getopt("n_dist") == "uniform")
        n = Grid(pms->getval("Nx"),pms->getval("Ny"),n_max);
    else assert(false && "Density distribution options are: <gaussian>, <exponential>, or <uniform>.");
    // option for adding a central ion hole to the x axis of the plasma
    if (pms->getopt("n_hole") == "true"){
        double hole_amp = pms->getval("n_hole_amp");
        double hole_min = 0;
        double hole_size = pms->getval("n_hole_size");
        Grid n_hole = Grid::Gaussian2D(getgrid("pos_x"),getgrid("pos_y"),hole_amp,hole_min,hole_size,pms->getval("sig_y"),0,0);
        n -= n_hole;
    }
    if (pms->getopt("n_iaw") == "true"){
        double iaw_amp = pms->getval("n_iaw_amp");
        double iaw_sig = pms->getval("n_iaw_sig");
        double k = 2*PI/iaw_sig;
        Grid x = getgrid("pos_x");
        Grid n_iaw = iaw_amp*n*(k*x).sin();
        n -= n_iaw;
    }
    return n;
}