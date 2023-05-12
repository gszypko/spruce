#include "module.hpp"
#include "plasmadomain.hpp"
#include "radiativelosses.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <cmath>
#include <iostream>
#include <memory>
#include <cassert>

RadiativeLosses::RadiativeLosses(PlasmaDomain &pd): Module(pd){
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    avg_losses = Grid::Zero(1,1);
}

void RadiativeLosses::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "cutoff_ramp") cutoff_ramp = std::stod(this_rhs);
        else if(this_lhs == "cutoff_temp") cutoff_temp = std::stod(this_rhs);
        else if(this_lhs == "epsilon") epsilon = std::stod(this_rhs);
        else if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "time_integrator") time_integrator = this_rhs;
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void RadiativeLosses::setupModule(){
    if(time_integrator == "") time_integrator == "euler";
    assert((time_integrator == "euler" || time_integrator == "rk2" || time_integrator == "rk4")
            && "Invalid time integrator given for Thermal Conduction module");
    if(output_to_file) avg_losses = Grid(m_pd.m_xdim,m_pd.m_ydim,0.0);
}

void RadiativeLosses::preIterateModule(double dt){
    curr_num_subcycles = numberSubcycles(dt);
}

void RadiativeLosses::iterateModule(double dt){
    Grid thermal_energy = m_pd.m_eqs->grid(IdealMHD::thermal_energy);
    Grid temp = m_pd.m_eqs->grid(IdealMHD::temp);
    Grid n = m_pd.m_eqs->grid(IdealMHD::n);
    Grid old_thermal_energy;
    if(output_to_file) old_thermal_energy = thermal_energy;
    double dt_subcycle = (dt/(double)curr_num_subcycles);

    if(time_integrator == "euler") for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
        thermal_energy = thermal_energy - m_pd.m_ghost_zone_mask*dt_subcycle*computeLosses(temp,n);
        thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
        temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*n)).max(m_pd.temp_min);
    }
    else if(time_integrator == "rk2") {
        Grid half_step, half_step_temp;
        for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
            half_step = thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*computeLosses(temp,n);
            half_step = half_step.max(m_pd.thermal_energy_min);
            half_step_temp = ((m_pd.m_adiabatic_index - 1.0)*half_step/(2*K_B*n)).max(m_pd.temp_min);
            
            thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*computeLosses(half_step_temp,n);
            thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
            temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*n)).max(m_pd.temp_min);
        }
    }
    else if(time_integrator == "rk4") {
        Grid intermediate, intermediate_temp, k1, k2, k3, k4;
        for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
            k1 = computeLosses(temp,n);

            intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*k1).max(m_pd.thermal_energy_min);
            intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*n)).max(m_pd.temp_min);
            k2 = computeLosses(intermediate_temp,n);

            intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*k2).max(m_pd.thermal_energy_min);
            intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*n)).max(m_pd.temp_min);
            k3 = computeLosses(intermediate_temp,n);

            intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*k3).max(m_pd.thermal_energy_min);
            intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*n)).max(m_pd.temp_min);
            k4 = computeLosses(intermediate_temp,n);

            thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
            thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
            temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*n)).max(m_pd.temp_min);
        }
    }


    if(output_to_file) avg_losses = (thermal_energy - old_thermal_energy)/dt;
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) = thermal_energy;
    m_pd.m_eqs->propagateChanges();
}

Grid RadiativeLosses::computeLosses() const {
    return computeLosses(m_pd.m_eqs->grid(IdealMHD::temp), m_pd.m_eqs->grid(IdealMHD::n));
}

//Compute the volumetric rate of energy loss to radiation
//using piecewise approximation for optically thin corona
//Result is stored in next_losses and is positive, such that change in energy is -1*next_losses
Grid RadiativeLosses::computeLosses(const Grid &temp, const Grid &n) const {
    Grid result = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    #pragma omp parallel for collapse(2)
    for (int i = m_pd.m_xl; i <= m_pd.m_xu; i++){
        for(int j = m_pd.m_yl; j <= m_pd.m_yu; j++){
            if(temp(i,j) < cutoff_temp) result(i,j) = 0.0;
            else {
                double logtemp = std::log10(temp(i,j));
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
                result(i,j) = std::pow(n(i,j),2.0)*chi*std::pow(temp(i,j),alpha);
                if(temp(i,j) < cutoff_temp + cutoff_ramp){
                    double ramp = (temp(i,j) - cutoff_temp)/cutoff_ramp;
                    result(i,j) *= ramp;
                }
            }
        }
    }
    return result;
}

//Computes number of subcycles necessary for the current iteration of the module
int RadiativeLosses::numberSubcycles(double dt){
    Grid losses = computeLosses();
    if(losses.max() == 0.0) return 0;
    double subcycle_dt = epsilon*(m_pd.m_eqs->grid(IdealMHD::thermal_energy)/losses).abs().min();
    return (int)(dt/subcycle_dt) + 1;
}


std::string RadiativeLosses::commandLineMessage() const
{
    return "Radiative Subcycles: " + std::to_string(curr_num_subcycles);
}

void RadiativeLosses::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if (output_to_file) {
        if (avg_losses.size() == 1) avg_losses = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
        var_names.push_back("rad");
        var_grids.push_back(avg_losses);
    }
}
