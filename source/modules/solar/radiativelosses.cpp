//radiativelosses.hpp
//Header for the Radiative Losses Module,
//an implementation of the abstract Module class
//Applies a piecewise power-law approximation
//of optically-thin radiation in the solar corona

#include "module.hpp"
#include "plasmadomain.hpp"
#include "radiativelosses.hpp"

RadiativeLosses::RadiativeLosses(PlasmaDomain &pd): Module(pd) {}

void RadiativeLosses::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "cutoff_ramp") cutoff_ramp = std::stod(this_rhs);
        else if(this_lhs == "cutoff_temp") cutoff_temp = std::stod(this_rhs);
        else if(this_lhs == "epsilon") epsilon = std::stod(this_rhs);
        else if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void RadiativeLosses::preIterateModule(double dt){
    next_losses = Grid(m_pd.m_xdim,m_pd.m_ydim,0.0);
    if(output_to_file) avg_losses = Grid(m_pd.m_xdim,m_pd.m_ydim,0.0);
    curr_num_subcycles = numberSubcycles(dt);
}

void RadiativeLosses::iterateModule(double dt){
    for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
        computeLosses();
        if(output_to_file) avg_losses += next_losses/curr_num_subcycles;
        m_pd.m_grids[PlasmaDomain::thermal_energy] = m_pd.m_grids[PlasmaDomain::thermal_energy] - m_pd.m_ghost_zone_mask*(dt/(double)curr_num_subcycles)*next_losses;
        m_pd.propagateChanges();
    }
}

//Compute the volumetric rate of energy loss to radiation
//using piecewise approximation for optically thin corona
//Result is stored in next_losses and is positive, such that change in energy is -1*next_losses
void RadiativeLosses::computeLosses(){
    assert(next_losses.rows() == m_pd.m_xdim && next_losses.cols() == m_pd.m_ydim && "This function assumes that the next_losses Grid has already been allocated");
    #pragma omp parallel for collapse(2)
    for (int i = m_pd.m_xl; i <= m_pd.m_xu; i++){
        for(int j = m_pd.m_yl; j <= m_pd.m_yu; j++){
            if(m_pd.m_grids[PlasmaDomain::temp](i,j) < cutoff_temp) next_losses(i,j) = 0.0;
            else {
                double logtemp = std::log10(m_pd.m_grids[PlasmaDomain::temp](i,j));
                double n = m_pd.m_grids[PlasmaDomain::rho](i,j)/m_pd.m_ion_mass;
                double chi, alpha;
                if(logtemp <= 4.97){
                    chi = 1.09e-31; //also adjust chi to ensure continuity
                    alpha = 2.0; //alpha 3 might be better approx?
                    // chi = 1.17e-36;
                    // alpha = 3.0;
                } else if(logtemp <= 5.67){
                    chi = 8.87e-17;
                    alpha = -1.0;
                } else if(logtemp <= 6.18){
                    chi = 1.90e-22;
                    alpha = 0.0;
                } else if(logtemp <= 6.55){
                    chi = 3.53e-13;
                    alpha = -1.5;
                } else if(logtemp <= 6.90){
                    chi = 3.46e-25;
                    alpha = 1.0/3.0;
                } else if(logtemp <= 7.63){
                    chi = 5.49e-16;
                    alpha = -1.0;
                } else{
                    chi = 1.96e-27;
                    alpha = 0.5;
                }
                next_losses(i,j) = n*n*chi*std::pow(m_pd.m_grids[PlasmaDomain::temp](i,j),alpha);
                if(m_pd.m_grids[PlasmaDomain::temp](i,j) < cutoff_temp + cutoff_ramp){
                    double ramp = (m_pd.m_grids[PlasmaDomain::temp](i,j) - cutoff_temp)/cutoff_ramp;
                    next_losses(i,j) *= ramp;
                }
            }
        }
    }
}

//Computes number of subcycles necessary for the current iteration of the module
int RadiativeLosses::numberSubcycles(double dt){
    computeLosses();
    if(next_losses.max() == 0.0) return 0;
    double subcycle_dt = epsilon*(m_pd.m_grids[PlasmaDomain::thermal_energy]/next_losses).abs().min();
    return (int)(dt/subcycle_dt) + 1;
}


std::string RadiativeLosses::commandLineMessage() const
{
    return "Radiative Subcycles: " + std::to_string(curr_num_subcycles);
}

void RadiativeLosses::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) const
{
    var_names.push_back("rad");
    var_grids.push_back(avg_losses);
}
