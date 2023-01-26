#include "viscosity.hpp"

Viscosity::Viscosity(PlasmaDomain &pd): Module(pd) {}

void Viscosity::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "visc_output_visc") m_output_visc = (rhs[i] == "true");
        else if(lhs[i] == "visc_output_lap") m_output_lap = (rhs[i] == "true");
        else if(lhs[i] == "visc_output_strength") m_output_strength = (rhs[i] == "true");
        else if(lhs[i] == "visc_output_timescale") m_output_timescale = (rhs[i] == "true");
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
    m_grids_lap = std::vector<Grid>(m_num_terms,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    m_grids_strength = std::vector<Grid>(m_num_terms,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    m_grids_dqdt = std::vector<Grid>(m_num_terms,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    m_dt_names.resize(m_num_terms);
    m_lap_names.resize(m_num_terms);
    m_strength_names.resize(m_num_terms);
    m_dqdt_names.resize(m_num_terms);
    for (int i=0; i<m_num_terms; i++){
        m_dt_names[i] = m_vars_to_evol[i] + "_dt";
        m_lap_names[i] = m_vars_to_evol[i] + "_lap";
        m_strength_names[i] = m_vars_to_evol[i] + "_str";
        m_dqdt_names[i] = m_vars_to_evol[i] + "_dqdt";
    } 
}

void Viscosity::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt)
{
    constructViscosityGrids(grids);
    for (int i=0; i<m_num_terms; i++)
        grids_dt[m_pd.m_eqs->name2evolvedindex(m_vars_to_evol[i])] += m_grids_dqdt[i]*m_pd.m_ghost_zone_mask;
        
}

void Viscosity::constructViscosityGrids(const std::vector<Grid>& grids)
{
    for (int i=0; i<m_num_terms; i++){
        // reference to plasma domain grids
        const Grid& dx = m_pd.grid("d_x");
        const Grid& dy = m_pd.grid("d_y");

        // determine viscosity strength factor
        if (m_visc_opt[i] == "local" || m_visc_opt[i] == "global") m_grids_strength[i] = Grid(m_pd.m_xdim,m_pd.m_ydim,m_strength[i]);
        else if (m_visc_opt[i] == "boundary") m_grids_strength[i] = getBoundaryViscosity(m_strength[i],m_length[i]);
        else assert("Viscosity option must be global, local, or boundary.");

        // determine timescale to base viscosity off of
        bool is_1F_model = m_pd.m_eqs->timescale().size() == 1;
        Grid dt; // grid that holds evolution timescale for species that viscosity term is being applied to
        if (is_1F_model) dt = m_pd.m_eqs->grid("dt"); // the main "dt" grid is used if 1F model
        else{ // if not 1F model, use timescale() indices to determine electron/ion timescale, as necessary
            if (m_species[i] == "i") dt = m_pd.m_eqs->grid(m_pd.m_eqs->timescale()[0]);
            else if (m_species[i] == "e") dt = m_pd.m_eqs->grid(m_pd.m_eqs->timescale()[1]);
            else assert(false && "Species must be <i> or <e>.");
        }
        double dt_min = dt.min(m_pd.m_xl_dt,m_pd.m_yl_dt,m_pd.m_xu_dt,m_pd.m_yu_dt);
        if (m_visc_opt[i] == "boundary" || m_visc_opt[i] == "local") m_grids_dt[i] = dt;
        else if (m_visc_opt[i] == "global") m_grids_dt[i] = Grid(m_pd.m_xdim,m_pd.m_ydim,dt_min);
        else assert("Viscosity option must be global, local, or boundary.");

        // compute viscosity coefficient and impose maximum
        Grid visc_coeff = m_grids_strength[i]*(dx.square()+dy.square())/2./m_grids_dt[i];
        double rk_timestep = m_pd.epsilon*m_pd.m_eqs->getDT().min(m_pd.m_xl_dt,m_pd.m_yl_dt,m_pd.m_xu_dt,m_pd.m_yu_dt);
        Grid visc_max = (dx.square()+dy.square())/2./rk_timestep;
        visc_coeff.min(visc_max);

        // compute Laplacian of variable to differentiate
        const Grid& grid_to_diff = grids[m_pd.m_eqs->name2index(m_vars_to_diff[i])];
        m_grids_lap[i] = m_pd.laplacian(grid_to_diff);
        
        // apply scaling of viscosity grid when applying viscosity to momentum density
        Grid scale_fac = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);
        if (is_momentum(m_vars_to_evol[i])){
            if (is_velocity(m_vars_to_diff[i])){
                // handle scaling for ions
                if (m_species[i] == "i"){
                    // ion scaling depends on whether 1F or 2F model
                    if (m_pd.m_eqs->is_var("i_n")) scale_fac = grids[m_pd.m_eqs->name2index("i_n")]*m_pd.m_ion_mass;
                    else scale_fac = grids[m_pd.m_eqs->name2index("n")]*m_pd.m_ion_mass;
                } // handle scaling for electrons
                else if (m_species[i] == "e") scale_fac = grids[m_pd.m_eqs->name2index("e_n")]*M_ELECTRON;
                else assert("Species must be e or i when applying viscosity to momentum.");
            }
            else{
                if(!is_momentum(m_vars_to_diff[i])) std::cerr << "Variable being differentiated is not consistent with variable viscosity is applied to.\n";
            }
        }

        // compute final viscosity term
        m_grids_dqdt[i] = visc_coeff*m_grids_lap[i]*scale_fac;

    }
}

Grid Viscosity::getBoundaryViscosity(double strength,double length) const
{
    // initialize grids and references to PlasmaDomain grids
    Grid result = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    const Grid& x = m_pd.grid("pos_x");
    const Grid& y = m_pd.grid("pos_y");
    // left boundary
    result += (-2.3*((x - x.min())/length).square()).exp()*strength;
    // right boundary
    result += (-2.3*((x - x.max())/length).square()).exp()*strength;
    // top boundary
    result += (-2.3*((y - y.max())/length).square()).exp()*strength;
    // bottom boundary
    result += (-2.3*((y - y.min())/length).square()).exp()*strength;
    // call to function min ensures that viscosity value cannot exceed strength
    // - this call is necessary because the boundary viscosity term is formed via superposition of many terms
    return result.min(strength);
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
    if (m_output_visc) {
        for (int i = 0; i<m_num_terms; i++){
            var_names.push_back(m_dqdt_names[i]);
            var_grids.push_back(m_grids_dqdt[i]);
        }
    }
    if (m_output_lap){
        for (int i = 0; i<m_num_terms; i++){
            var_names.push_back(m_lap_names[i]);
            var_grids.push_back(m_grids_lap[i]);
        }
    }
    if (m_output_strength){
        for (int i = 0; i<m_num_terms; i++){
            var_names.push_back(m_strength_names[i]);
            var_grids.push_back(m_grids_strength[i]);
        }
    }
    if (m_output_timescale){
        for (int i = 0; i<m_num_terms; i++){
            var_names.push_back(m_dt_names[i]);
            var_grids.push_back(m_grids_dt[i]);
        }
    }  
}