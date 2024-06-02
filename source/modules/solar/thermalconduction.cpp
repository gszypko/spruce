#include "module.hpp"
#include "plasmadomain.hpp"
#include "thermalconduction.hpp"
#include "constants.hpp"
#include "idealmhd.hpp"
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>

ThermalConduction::ThermalConduction(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
}

void ThermalConduction::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "flux_saturation")  flux_saturation = (this_rhs == "true");
        else if(this_lhs == "epsilon")  epsilon = std::stod(this_rhs);
        else if(this_lhs == "dt_subcycle_min")  dt_subcycle_min = std::stod(this_rhs);
        else if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "time_integrator") time_integrator = this_rhs;
        else if(this_lhs == "inactive_mode") inactive_mode = (this_rhs == "true");
        else if(this_lhs == "weakening_factor")  weakening_factor = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized for Thermal Conduction Module.\n";
    }
}

void ThermalConduction::setupModule(){
    if(output_to_file) avg_change = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    else avg_change = Grid::Zero(1,1);
    if(time_integrator == "") time_integrator = "euler";
    assert((time_integrator == "euler" || time_integrator == "rk2" || time_integrator == "rk4")
            && "Invalid time integrator given for Thermal Conduction module");
}

void ThermalConduction::preIterateModule(double dt){
    curr_num_subcycles = numberSubcycles(dt);
}

void ThermalConduction::iterateModule(double dt){
    //Subcycle to simulate field-aligned thermal conduction
    Grid thermal_energy = m_pd.m_eqs->grid(IdealMHD::thermal_energy), temp = m_pd.m_eqs->grid(IdealMHD::temp), m_n = m_pd.m_eqs->grid(IdealMHD::n);
    std::vector<Grid> b_hat = {m_pd.m_eqs->grid(IdealMHD::b_hat_x), m_pd.m_eqs->grid(IdealMHD::b_hat_y)};
    // Grid thermal_energy_next;
    if(output_to_file) avg_change = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    double dt_subcycle = (dt/(double)curr_num_subcycles);
    Grid old_thermal_energy = thermal_energy;
    if(time_integrator == "euler") for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
        thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*thermalEnergyDerivative(temp,b_hat);
        thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
        temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2.0*K_B*m_n)).max(m_pd.temp_min);
    }
    else if(time_integrator == "rk2") {
        Grid half_step, half_step_temp;
        for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
            half_step = thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*thermalEnergyDerivative(temp,b_hat);
            half_step = half_step.max(m_pd.thermal_energy_min);
            half_step_temp = ((m_pd.m_adiabatic_index - 1.0)*half_step/(2*K_B*m_n)).max(m_pd.temp_min);
            
            thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*thermalEnergyDerivative(half_step_temp,b_hat);
            thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
            temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*m_n)).max(m_pd.temp_min);
        }
    }
    else if(time_integrator == "rk4") {
        Grid intermediate, intermediate_temp, k1, k2, k3, k4;
        for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
            k1 = thermalEnergyDerivative(temp,b_hat);

            intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*k1).max(m_pd.thermal_energy_min);
            intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*m_n)).max(m_pd.temp_min);
            k2 = thermalEnergyDerivative(intermediate_temp,b_hat);

            intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*k2).max(m_pd.thermal_energy_min);
            intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*m_n)).max(m_pd.temp_min);
            k3 = thermalEnergyDerivative(intermediate_temp,b_hat);

            intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*k3).max(m_pd.thermal_energy_min);
            intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*m_n)).max(m_pd.temp_min);
            k4 = thermalEnergyDerivative(intermediate_temp,b_hat);

            thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
            thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
            temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*m_n)).max(m_pd.temp_min);
        }
    }
    if(output_to_file) avg_change = (thermal_energy - old_thermal_energy)/dt;
    if(inactive_mode) return;
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) = thermal_energy;
    m_pd.m_eqs->propagateChanges();
}

Grid ThermalConduction::thermalEnergyDerivative(Grid &m_temp, std::vector<Grid> b_hat) const{
    Grid dtemp_dx = m_pd.derivative1D(m_temp,0);
    Grid dtemp_dy = m_pd.derivative1D(m_temp,1);
    std::vector<Grid> gradt = {dtemp_dx,dtemp_dy};
    Grid b_dot_gradt = Grid::DotProduct2D(b_hat,gradt);
    std::vector<Grid> term1 = {b_hat[0]*m_pd.secondDerivative1D(m_temp,0) + b_hat[1]*m_pd.derivative1D(dtemp_dx,1),
                                b_hat[0]*m_pd.derivative1D(dtemp_dy,0) + b_hat[1]*m_pd.secondDerivative1D(m_temp,1)};
    std::vector<Grid> term2 = {dtemp_dx*m_pd.derivative1D(b_hat[0],0) + dtemp_dy*m_pd.derivative1D(b_hat[0],1),
                                dtemp_dx*m_pd.derivative1D(b_hat[1],0) + dtemp_dy*m_pd.derivative1D(b_hat[1],1)};
    std::vector<Grid> term3 = Grid::CrossProduct2DZ(gradt,m_pd.curl2D(b_hat[0],b_hat[1]));
    std::vector<Grid> term_total(2,Grid::Zero(1,1));
    for(int i : {0,1}) term_total[i] = 5.0/2.0*m_temp.pow(3.0/2.0)*gradt[i]*b_dot_gradt + m_temp.pow(5.0/2.0)*(term1[i] + term2[i] + term3[i]);
    Grid result = (weakening_factor*KAPPA_0)*((m_temp.pow(5.0/2.0) * b_dot_gradt * m_pd.divergence2D(b_hat)) + Grid::DotProduct2D(b_hat,term_total));
    if(flux_saturation) {
        std::vector<Grid> sat_terms = saturationTerms();
        result = sat_terms[0]*result + sat_terms[1];
    }
    return result;
}

//Computes number of subcycles necessary for the current iteration of the module
int ThermalConduction::numberSubcycles(double dt){
    Grid dt_subcycle;
    if(!flux_saturation){
        // This time-step calculation comes from (nkT)/dT = (2/7 kappa T^7/2)/dx^2
        dt_subcycle = K_B/(weakening_factor*KAPPA_0)*(m_pd.m_eqs->grid(IdealMHD::rho)/m_pd.m_ion_mass)*m_pd.m_grids[PlasmaDomain::d_x]*m_pd.m_grids[PlasmaDomain::d_y]/m_pd.m_eqs->grid(IdealMHD::temp).pow(2.5);
    } else {
        Grid field_temp_gradient = m_pd.derivative1D(m_pd.m_eqs->grid(IdealMHD::temp),0)*m_pd.m_eqs->grid(IdealMHD::b_hat_x)
                                + m_pd.derivative1D(m_pd.m_eqs->grid(IdealMHD::temp),1)*m_pd.m_eqs->grid(IdealMHD::b_hat_y);
        if(field_temp_gradient.abs().max() == 0.0) return 0;
        Grid kappa_modified = saturatedKappa();
        dt_subcycle = K_B/kappa_modified*(m_pd.m_eqs->grid(IdealMHD::rho)/m_pd.m_ion_mass)*m_pd.m_grids[PlasmaDomain::d_x]*m_pd.m_grids[PlasmaDomain::d_y];
    }
    double min_dt_subcycle = std::max(epsilon*dt_subcycle.min(m_pd.m_xl,m_pd.m_yl,m_pd.m_xu,m_pd.m_yu),dt_subcycle_min);
    return (int)(dt/min_dt_subcycle) + 1;
}

//Computes 1D cell-centered conductive flux from temperature "temp"
//Flux computed in direction indicated by "index": 0 for x, 1 for y
//k0 is conductive coefficient
Grid ThermalConduction::oneDimConductiveFlux(const Grid &temp, const Grid &rho, double k0, int index) const {
    Grid kappa_max = m_pd.m_grids[PlasmaDomain::d_x]*m_pd.m_grids[PlasmaDomain::d_y]*K_B*(rho/m_pd.m_ion_mass)/dt_subcycle_min;
    return -(k0*temp.pow(5.0/2.0)).min(kappa_max)*m_pd.derivative1D(temp,index);
}

//Computes cell-centered, field-aligned conductive flux from temperature "temp"
//temp is temperature Grid
//b_hat_x, b_hat_y are the components of the *unit* vector b_hat
//k0 is conductive coefficient
//Output is written to flux_out_x and flux_out_y
void ThermalConduction::fieldAlignedConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &temp, const Grid &rho,
                                    const Grid &b_hat_x, const Grid &b_hat_y, double k0) const {
    int xdim = temp.rows();
    int ydim = temp.cols();
    Grid con_flux_x = oneDimConductiveFlux(temp, rho, k0, 0);
    Grid con_flux_y = oneDimConductiveFlux(temp, rho, k0, 1);
    #pragma omp parallel for collapse(2)
    for (int i = m_pd.m_xl; i <= m_pd.m_xu; i++){
        for(int j = m_pd.m_yl; j <= m_pd.m_yu; j++){
            double flux_magnitude = con_flux_x(i,j)*b_hat_x(i,j) + con_flux_y(i,j)*b_hat_y(i,j);
            flux_out_x(i,j) = flux_magnitude*b_hat_x(i,j);
            flux_out_y(i,j) = flux_magnitude*b_hat_y(i,j);
        }
    }
}

//Computes saturated conductive flux at each point in grid,
//then ensures that provided fluxes do not exceed the saturation point
void ThermalConduction::saturateConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &rho, const Grid &temp) const {
    Grid sat_flux_mag = (1.0/6.0)*(3.0/2.0)*(rho/m_pd.m_ion_mass)*(K_B*temp).pow(1.5)/std::sqrt(M_ELECTRON);
    Grid flux_mag = (flux_out_x.square() + flux_out_y.square()).sqrt();
    Grid scale_factor = sat_flux_mag /((sat_flux_mag.square() + flux_mag.square()).sqrt());
    flux_out_x *= scale_factor;
    flux_out_y *= scale_factor;
}

//NOTE: This computes the temperature-dependent kappa factor, i.e. kappa_0*temp^5/2 in the unsaturated case
//such that the heat flux is kappa * grad(T)
Grid ThermalConduction::saturatedKappa() const {
    Grid kappa_modified(m_pd.m_xdim,m_pd.m_ydim,0.0);
    Grid con_flux_x(m_pd.m_xdim,m_pd.m_ydim,0.0);
    Grid con_flux_y(m_pd.m_xdim,m_pd.m_ydim,0.0);
    Grid field_temp_gradient = m_pd.derivative1D(m_pd.m_eqs->grid(IdealMHD::temp),0)*m_pd.m_eqs->grid(IdealMHD::b_hat_x)
                            + m_pd.derivative1D(m_pd.m_eqs->grid(IdealMHD::temp),1)*m_pd.m_eqs->grid(IdealMHD::b_hat_y);
    fieldAlignedConductiveFlux(con_flux_x, con_flux_y, m_pd.m_eqs->grid(IdealMHD::temp), m_pd.m_eqs->grid(IdealMHD::rho),
                                m_pd.m_eqs->grid(IdealMHD::b_hat_x), m_pd.m_eqs->grid(IdealMHD::b_hat_y), (weakening_factor*KAPPA_0));
    saturateConductiveFlux(con_flux_x, con_flux_y, m_pd.m_eqs->grid(IdealMHD::rho), m_pd.m_eqs->grid(IdealMHD::temp));
    Grid flux_mag = (con_flux_x.square() + con_flux_y.square()).sqrt();
    kappa_modified = (flux_mag/field_temp_gradient).abs();
    return kappa_modified;
}

//NOTE: This computes the multiplicative factor that is equivalent to applying flux saturation to the flux
//Returns a vector containing {saturationCoefficient,saturationAdditiveFactor}
//saturationCoefficient is the multiplicative factor to multiply into the computed energy ROC
//saturationAdditiveFactor is the quantity to add to the compute energy ROC
//Applying both of the above to the computed energy ROC is equivalent to applying flux saturation to the heat flux
std::vector<Grid> ThermalConduction::saturationTerms() const {
    Grid con_flux_x(m_pd.m_xdim,m_pd.m_ydim,0.0);
    Grid con_flux_y(m_pd.m_xdim,m_pd.m_ydim,0.0);
    fieldAlignedConductiveFlux(con_flux_x, con_flux_y, m_pd.m_eqs->grid(IdealMHD::temp), m_pd.m_eqs->grid(IdealMHD::rho),
                                m_pd.m_eqs->grid(IdealMHD::b_hat_x), m_pd.m_eqs->grid(IdealMHD::b_hat_y), (weakening_factor*KAPPA_0));
    Grid flux_mag = (con_flux_x.square() + con_flux_y.square()).sqrt();
    std::vector<Grid> raw_flux = {con_flux_x,con_flux_y};
    saturateConductiveFlux(con_flux_x, con_flux_y, m_pd.m_eqs->grid(IdealMHD::rho), m_pd.m_eqs->grid(IdealMHD::temp));
    Grid sat_flux_mag = (con_flux_x.square() + con_flux_y.square()).sqrt();
    Grid coefficient(m_pd.m_xdim,m_pd.m_ydim,0.0);
    for(int i=0; i<m_pd.m_xdim; i++) for (int j=0; j<m_pd.m_ydim; j++){
        if (flux_mag(i,j) != 0.0) coefficient(i,j) = sat_flux_mag(i,j)/flux_mag(i,j);
        else coefficient(i,j) = 1.0;
    }
    Grid additive_factor = -1.0*m_pd.m_ghost_zone_mask*Grid::DotProduct2D({m_pd.derivative1D(coefficient,0),m_pd.derivative1D(coefficient,1)},raw_flux);
    return {coefficient,additive_factor};
}

void ThermalConduction::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if (output_to_file) {
        if (avg_change.size() == 1) avg_change = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
        var_names.push_back("thermal_conduction");
        var_grids.push_back(avg_change);
    }
}


std::string ThermalConduction::commandLineMessage() const
{
    return "Thermal Subcycles: " + std::to_string(curr_num_subcycles) + (inactive_mode ? " (Not Applied)" : "");
}