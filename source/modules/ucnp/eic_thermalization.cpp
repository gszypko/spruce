#include "eic_thermalization.hpp"

#define LENGTH_DECAY [](double r,double l){return exp(-r/l);} 

EICThermalization::EICThermalization(PlasmaDomain &pd): Module(pd) {}

void EICThermalization::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "timescale") m_timescale = std::stod(rhs[i]);
        else if(lhs[i] == "lengthscale") m_lengthscale = std::stod(rhs[i]);
        else if(lhs[i] == "strength") m_strength = std::stod(rhs[i]);
        else if(lhs[i] == "output_to_file") output_to_file = (rhs[i] == "true");
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void EICThermalization::setupModule()
{
    // ensure that grid dependencies exist within EquationSet
    std::vector<std::string> grid_names = m_pd.m_eqs->def_var_names();
    for (auto grid : m_required_grids){
        auto it = std::find(grid_names.begin(),grid_names.end(),grid);
        if (it==grid_names.end()){
            std::cerr << "Grid <" << grid << "> was not found within the EquationSet." << std::endl;
            assert(false);
        }
        int loc = std::distance(grid_names.begin(),it);
        m_grid_ind[grid] = loc;
    }
    // initialize grids
    m_Fcx = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    m_Fcy = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    m_dPx = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    m_dPy = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
}

// returns rho_c values in cgs at each distance in <r>
Grid EICThermalization::compute_charge_density(const Grid& r, const Grid& n) const
{
    double time_decay = exp(-m_pd.m_time/m_timescale);
    Grid length_decay = r.for_each(m_lengthscale,[](double a,double b){return exp(-a/b);});
    Grid rho_c = (m_strength*time_decay*E)*n*length_decay;
    return rho_c;
}

Grid EICThermalization::compute_total_charge(const Grid& r,const Grid& rho_c) const
{
    return (4.*PI)*Grid::trapzcum(r,r.square()*rho_c);
}

void EICThermalization::postIterateModule(double dt)
{
    if (m_pd.m_time < 3*m_timescale){
    // references to 2D grids
    const Grid& x {m_pd.m_internal_grids[PlasmaDomain::pos_x]};
    const Grid& y {m_pd.m_internal_grids[PlasmaDomain::pos_y]};
    const Grid& n {m_pd.grid(m_grid_ind.at("n"))};
    const Grid& press {m_pd.grid(m_grid_ind.at("press"))};
    Grid& mom_x {m_pd.grid(m_grid_ind.at("mom_x"))};
    Grid& mom_y {m_pd.grid(m_grid_ind.at("mom_y"))};
    // compute distance from plasma center
    Grid r_sq = x.square()+y.square();
    Grid r = r_sq.sqrt();
    // bin the distances and densities
    Grid r_bin;
    Grid n_bin = Grid::bin_as_list(r,n,101,r_bin);
    // compute the charge density
    m_dPx = m_pd.derivative1D(press, 0);
    m_dPy = m_pd.derivative1D(press, 1);
    Grid rho_c = compute_charge_density(r_bin,n_bin);
    Grid Q_vec = compute_total_charge(r_bin,rho_c);
    // compute electric field
    Grid Q = Grid::interp_as_list(r_bin,Q_vec,r);
    Grid r_cubed = r_sq*r;
    Grid Ex = x*Q/r_cubed;
    Grid Ey = y*Q/r_cubed;
    // compute force profile
    m_Fcx = E*n*Ex;
    m_Fcy = E*n*Ey;
    // add forces
    mom_x += m_Fcx*dt;
    mom_y += m_Fcy*dt;
    m_pd.m_eqs->propagateChanges();
    }
}

std::string EICThermalization::commandLineMessage() const
{
    return "Coulomb Explosion On";
}

void EICThermalization::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if (output_to_file) {
        var_names.push_back("Fcx");
        var_grids.push_back(m_Fcx);
        var_names.push_back("Fcy");
        var_grids.push_back(m_Fcy);
        var_names.push_back("dPx");
        var_grids.push_back(m_dPx);
        var_names.push_back("dPy");
        var_grids.push_back(m_dPy);
    }
}