#include "viscosity.hpp"

Viscosity::Viscosity(PlasmaDomain &pd): Module(pd) {}

void Viscosity::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "visc_output_to_file") m_output_to_file = (rhs[i] == "true");
        else if(lhs[i] == "visc_opt") m_inp_visc_opt = rhs[i];
        else if(lhs[i] == "visc_strength") m_inp_strength = rhs[i];
        else if(lhs[i] == "visc_vars_to_diff") m_inp_vars_to_diff = rhs[i];
        else if(lhs[i] == "visc_vars_to_evol") m_inp_vars_to_evol = rhs[i];
        else if(lhs[i] == "visc_length") m_inp_length = rhs[i];
        else if(lhs[i] == "visc_species") m_inp_species = rhs[i];
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void Viscosity::setupModule()
{
    // clear white space from inputs
    clearWhitespace(m_inp_visc_opt);
    clearWhitespace(m_inp_strength);
    clearWhitespace(m_inp_vars_to_diff);
    clearWhitespace(m_inp_vars_to_evol);
    clearWhitespace(m_inp_length);

    // parse comma-delimited inputs
    std::vector<std::string> temp;
    temp = splitString(m_inp_visc_opt,',');
    m_num_terms = temp.size();

    if (m_num_terms == 0) assert("No inputs detected."); // abort process if input strings are empty

    m_visc_opt.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++) m_visc_opt[i] = temp[i];
    
    temp = splitString(m_inp_strength,',');
    assert(m_num_terms==temp.size());
    m_strength.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++) m_strength[i] = stod(temp[i]);
    
    temp = splitString(m_inp_vars_to_diff,',');
    assert(m_num_terms==temp.size());
    m_vars_to_diff.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++) m_vars_to_diff[i] = temp[i];

    temp = splitString(m_inp_vars_to_evol,',');
    assert(m_num_terms==temp.size());
    m_vars_to_evol.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++) m_vars_to_evol[i] = temp[i];

    temp = splitString(m_inp_length,',');
    assert(m_num_terms==temp.size());
    m_length.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++) m_length[i] = stod(temp[i]);

    temp = splitString(m_inp_species,',');
    assert(m_num_terms==temp.size());
    m_species.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++) m_species[i] = temp[i];

    // ensure that options are compatible with current equation set
    for (auto val : m_strength) assert(val<=1 && "Strength constants must be less than or equal to one.");
    for (auto name : m_vars_to_diff) assert(m_pd.m_eqs->is_var(name) && "Each variable to differentiate must be a valid variable within the chosen equation set.");
    for (auto name : m_vars_to_evol) assert(m_pd.m_eqs->is_evolved_var(name) && "Each variable to evolve must be a valid evolved variable within the chosen equation set.");
    for (auto val : m_length) assert(val>=0 && "Length constants must greater than or equal to zero.");

    // initialize viscosity grids
    m_grids_dt = std::vector<Grid>(m_num_terms,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    m_grid_names.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++) m_grid_names[i] = m_vars_to_evol[i] + "_dt";
}

void Viscosity::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt)
{
    for (int i=0; i<m_num_terms; i++){
        m_grids_dt[i] = constructViscosityGrid(m_visc_opt[i], m_strength[i],m_vars_to_diff[i],m_vars_to_evol[i],m_length[i],m_species[i],grids);
        grids_dt[m_pd.m_eqs->name2evolvedindex(m_vars_to_evol[i])] += m_grids_dt[i]*m_pd.m_ghost_zone_mask;
    }
}

Grid Viscosity::constructViscosityGrid(std::string opt,double strength, std::string var_to_diff, std::string var_to_evol, double length, std::string species, const std::vector<Grid>& grids) const
{
    // reference to plasma domain grids
    const Grid& dx = m_pd.grid("d_x");
    const Grid& dy = m_pd.grid("d_y");
    const Grid& dt = grids[m_pd.m_eqs->name2index("dt")];
    const Grid& grid_to_diff = grids[m_pd.m_eqs->name2index(var_to_diff)];
    double dt_min = dt.min(m_pd.m_xl,m_pd.m_yl,m_pd.m_xu,m_pd.m_yu);
    // construct initial portion of viscosity grid of form eps*(dx^2+dy^2)/2/dt*Lap(q)
    Grid result {(dx.square()+dy.square())/2.};
    // apply strength constant depending on viscosity type
    if (opt == "local" || opt == "global") result *= strength;
    else if (opt == "boundary") result *= getBoundaryViscosity(strength,length);
    else assert("Viscosity option must be global, local, or boundary.");
    // apply time information depending on viscosity type
    if (opt == "local") result /= dt;
    else if (opt == "boundary" || opt == "global") result /= dt_min;
    else assert("Viscosity option must be global, local, or boundary.");
    // apply Laplacian to variable
    result *= m_pd.laplacian(grid_to_diff);
    // apply scaling of viscosity grid, depending on relation of evolved variable to differentiated variable
    if (is_momentum(var_to_evol) && is_velocity(var_to_diff)){
        // handle scaling for ions
        if (species == "i"){
            // ion scaling depends on whether 1F or 2F model
            if (m_pd.m_eqs->is_var("i_n")) result *= grids[m_pd.m_eqs->name2index("i_n")]*m_pd.m_ion_mass;
            else result *= grids[m_pd.m_eqs->name2index("n")]*m_pd.m_ion_mass;
        } // handle scaling for electrons
        else if (species == "e") result *= grids[m_pd.m_eqs->name2index("e_n")]*M_ELECTRON;
        else assert("Species must be e or i when applying viscosity to momentum.");
    }

    // // viscous forces
    // const Grid& d_x = m_pd.m_grids[PlasmaDomain::d_x];
    // const Grid& d_y = m_pd.m_grids[PlasmaDomain::d_y];
    // Grid global_visc_coeff = .1*0.5*(d_x.square()+d_y.square())/dt;
    // Grid result = global_visc_coeff*m_pd.laplacian(grid_to_diff);


    return result;
}

Grid Viscosity::getBoundaryViscosity(double strength,double length) const
{
    return Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
}

bool Viscosity::is_momentum(std::string name) const
{
    auto it = std::find(m_momenta.begin(),m_momenta.end(),name);
    return it != m_momenta.end();
}

bool Viscosity::is_velocity(std::string name) const
{
    auto it = std::find(m_velocities.begin(),m_velocities.end(),name);
    return it != m_velocities.end();
}

void Viscosity::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if (m_output_to_file) {
        for (int i = 0; i<m_num_terms; i++){
            var_names.push_back(m_grid_names[i]);
            var_grids.push_back(m_grids_dt[i]);
        }
    }
}