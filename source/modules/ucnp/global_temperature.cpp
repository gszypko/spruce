#include "global_temperature.hpp"

GlobalTemperature::GlobalTemperature(PlasmaDomain &pd): Module(pd) {}

void GlobalTemperature::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs)
{
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "gt_species") m_species = splitString(rhs[i],',');
        else if(lhs[i] == "gt_strength") m_strength = stod(rhs[i]);
        else if(lhs[i] == "gt_use_diffusion") m_use_diffusion = rhs[i]=="true";
        else if(lhs[i] == "gt_use_global_temp") m_use_global_temp = rhs[i]=="true";
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void GlobalTemperature::setupModule()
{   
    // process config inputs
    m_species_ind.resize(m_species.size(),-1);
    m_species_mass.resize(m_species.size());

    // check that requested species are valid
    for (int i=0; i<m_species.size(); i++){
        for (int j=0; j<m_pd.m_eqs->species().size(); j++){
            if (m_species[i] == m_pd.m_eqs->species()[j]) m_species_ind[i] = j;
        }
        if (m_species_ind[i] < 0){
            std::cerr << "Species <" << m_species[i] << "> does not correspond to a species within the active equation set." << std::endl;
            assert(false && "<gt_species> was specified incorrectly in the .config file.");
        }
    }

    // get mass for each species
    for (int i=0; i<m_species.size(); i++){
        if (m_species[i] == "e") m_species_mass[i] = M_ELECTRON;
        else if (m_species[i] == "i") m_species_mass[i] = m_pd.m_ion_mass;
        else assert(false && "Species name must be <e> or <i>.");
    }

    // compute spatial grid for diffusion term
    m_dr = (m_pd.grid("d_x").square()+m_pd.grid("d_y").square())*m_pd.m_ghost_zone_mask/2.;
}

void GlobalTemperature::preRecomputeDerivedModule(std::vector<Grid>& grids) const
{
    if (m_use_global_temp){
        // get vectors (i.e., 1D grid) for position along each axis
        Grid x_vec = m_pd.grid("pos_x").col(1);
        Grid y_vec = m_pd.grid("pos_y").row(1);
        for (int i=0; i<m_species.size(); i++){
            // grid references to needed equation set grids for current species
            const Grid& rho = grids[m_pd.m_eqs->densities()[m_species_ind[i]]];
            Grid& thermal_energy = grids[m_pd.m_eqs->thermal_energies()[m_species_ind[i]]];

            // compute global temperature
            double rho_int = Grid::Trapz2D(x_vec,y_vec,rho);
            double T_global = Grid::Trapz2D(x_vec,y_vec,thermal_energy)/rho_int; // this is effectively integrating rho*temp because thermal_energy \propto rho*temp
            
            // update grids to be consistent with global temperature
            thermal_energy = (T_global*rho).max(m_pd.thermal_energy_min); // scale global temperature back to thermal energy
        }
    }

}

void GlobalTemperature::postIterateModule(double dt)
{
    if (m_use_diffusion){
        for (int i=0; i<m_species.size(); i++){
            // determine diffusion timescale from strength factor
            double dt_fluid = dt/m_pd.epsilon;
            m_dt_diff = dt_fluid/m_strength;
            m_coeff = m_dr/m_dt_diff;
            m_dt_rk = m_dt_diff*m_pd.epsilon;
            m_num_rk_steps = dt/m_dt_rk; 

            // preallocate arrays
            Grid diff_term = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
            Grid temp_mid = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);

            // get equation set grids
            Grid& temp = m_pd.m_eqs->grid(m_pd.m_eqs->temperatures()[m_species_ind[i]]); 
            const Grid& n = m_pd.m_eqs->grid(m_pd.m_eqs->number_densities()[m_species_ind[i]]);
            Grid& thermal_energy = m_pd.m_eqs->grid(m_pd.m_eqs->thermal_energies()[m_species_ind[i]]);       
            for (int j=0; j<m_num_rk_steps; j++){
                diff_term = m_coeff*m_pd.laplacian(temp);
                temp_mid = temp + diff_term*m_dt_rk/2.;
                diff_term = m_coeff*m_pd.laplacian(temp_mid);
                temp += diff_term*m_dt_rk;
            }
            thermal_energy = n*K_B*temp/(m_pd.m_adiabatic_index-1.);
            m_pd.m_eqs->propagateChanges();        
        }
    }
}