#include "module.hpp"
#include "plasmadomain.hpp"
#include "anomalousresistivity.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include "solarutils.hpp"
#include <sstream>
#include <iostream>
#include <cmath>

AnomalousResistivity::AnomalousResistivity(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    time_scale = -1.0;
}

void AnomalousResistivity::setupModule(){
    assert(time_scale > 0.0 && "Anomalous Resistivity time_scale must be specified, and positive");
    diffusivity = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    joule_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    anomalous_template = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    if(metric_smoothing){
        assert(smoothing_sigma > 0.0 && "smoothing_sigma must be a positive number");
        kernel_radius = nearbyint(4.0*smoothing_sigma);
        int kernel_size = 2*kernel_radius + 1;
        smoothing_kernel = SolarUtils::GaussianGrid(kernel_size,kernel_size,0.0,1.0/(2.0*PI*smoothing_sigma*smoothing_sigma),
                                                    smoothing_sigma,smoothing_sigma,kernel_radius,kernel_radius);
    }
    else smoothing_kernel = Grid::Zero(1,1);
}

void AnomalousResistivity::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "time_scale") time_scale = std::stod(rhs[i]);
        else if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "metric_coeff") metric_coeff = std::stod(this_rhs);
        else if(this_lhs == "smoothing_sigma") smoothing_sigma = std::stod(this_rhs);
        else if(this_lhs == "safety_factor") safety_factor = std::stod(this_rhs);
        else if(this_lhs == "metric_smoothing") metric_smoothing = (this_rhs == "true");
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void AnomalousResistivity::iterateModule(double dt){
    computeDiffusion();
    Grid bi_x = m_pd.m_eqs->grid(IdealMHD::bi_x), bi_y = m_pd.m_eqs->grid(IdealMHD::bi_y), bi_z = m_pd.m_eqs->grid(IdealMHD::bi_z);
    Grid thermal_energy = m_pd.m_eqs->grid(IdealMHD::thermal_energy);
    // Grid thermal_energy_next;
    if(output_to_file) avg_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    double dt_subcycle = (dt/(double)curr_num_subcycles);
    Grid old_thermal_energy = thermal_energy;
    // if(time_integrator == "euler") for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
    for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
        // thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*thermalEnergyDerivative(temp,b_hat);
        // thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
        // temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2.0*K_B*m_n)).max(m_pd.temp_min);
        computeTemplate(m_pd.m_grids[PlasmaDomain::be_x]+bi_x,m_pd.m_grids[PlasmaDomain::be_y]+bi_y);
        Grid coeff = m_pd.m_ghost_zone_mask*anomalous_template*diffusivity;
        bi_x += dt_subcycle*coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_x])+m_pd.laplacian(bi_x));
        bi_y += dt_subcycle*coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_y])+m_pd.laplacian(bi_y));
        bi_z += dt_subcycle*coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_z])+m_pd.laplacian(bi_z));
        Grid electrical_resistivity = 4.0*PI/C/C*coeff;
        Grid current_density = C/(4.0*PI)*(m_pd.curl2D(m_pd.m_grids[PlasmaDomain::be_x]+bi_x,m_pd.m_grids[PlasmaDomain::be_y]+bi_y)).abs();
        joule_heating = electrical_resistivity*current_density*current_density;
        thermal_energy += dt_subcycle*joule_heating;
    }
    // else if(time_integrator == "rk2") {
    //     Grid half_step, half_step_temp;
    //     for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
    //         half_step = thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*thermalEnergyDerivative(temp,b_hat);
    //         half_step = half_step.max(m_pd.thermal_energy_min);
    //         half_step_temp = ((m_pd.m_adiabatic_index - 1.0)*half_step/(2*K_B*m_n)).max(m_pd.temp_min);
            
    //         thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*thermalEnergyDerivative(half_step_temp,b_hat);
    //         thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
    //         temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*m_n)).max(m_pd.temp_min);
    //     }
    // }
    // else if(time_integrator == "rk4") {
    //     Grid intermediate, intermediate_temp, k1, k2, k3, k4;
    //     for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
    //         k1 = thermalEnergyDerivative(temp,b_hat);

    //         intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*k1).max(m_pd.thermal_energy_min);
    //         intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*m_n)).max(m_pd.temp_min);
    //         k2 = thermalEnergyDerivative(intermediate_temp,b_hat);

    //         intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(0.5*dt_subcycle)*k2).max(m_pd.thermal_energy_min);
    //         intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*m_n)).max(m_pd.temp_min);
    //         k3 = thermalEnergyDerivative(intermediate_temp,b_hat);

    //         intermediate = (thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*k3).max(m_pd.thermal_energy_min);
    //         intermediate_temp = ((m_pd.m_adiabatic_index - 1.0)*intermediate/(2*K_B*m_n)).max(m_pd.temp_min);
    //         k4 = thermalEnergyDerivative(intermediate_temp,b_hat);

    //         thermal_energy = thermal_energy + m_pd.m_ghost_zone_mask*(dt_subcycle)*(k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;
    //         thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
    //         temp = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*m_n)).max(m_pd.temp_min);
    //     }
    // }
    if(output_to_file) avg_heating = (thermal_energy - old_thermal_energy)/dt;
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) = thermal_energy;
    m_pd.m_eqs->grid(IdealMHD::bi_x) = bi_x;
    m_pd.m_eqs->grid(IdealMHD::bi_y) = bi_y;
    m_pd.m_eqs->grid(IdealMHD::bi_z) = bi_z;
    m_pd.m_eqs->propagateChanges();
}


//std::vector<Grid> grids_dt {d_rho_dt,d_mom_x_dt,d_mom_y_dt,d_mom_z_dt,d_thermal_energy_dt,d_bi_x_dt,d_bi_y_dt,d_bi_z_dt};
void AnomalousResistivity::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt){
    computeDiffusion();
    computeTemplate(m_pd.m_grids[PlasmaDomain::be_x]+grids[IdealMHD::bi_x],m_pd.m_grids[PlasmaDomain::be_y]+grids[IdealMHD::bi_y]);
    Grid coeff = m_pd.m_ghost_zone_mask*anomalous_template*diffusivity;
    grids_dt[5] += coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_x])+m_pd.laplacian(grids[IdealMHD::bi_x]));
    grids_dt[6] += coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_y])+m_pd.laplacian(grids[IdealMHD::bi_y]));
    grids_dt[7] += coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_z])+m_pd.laplacian(grids[IdealMHD::bi_z]));
    Grid electrical_resistivity = 4.0*PI/C/C*coeff;
    Grid current_density = C/(4.0*PI)*(m_pd.curl2D(m_pd.m_grids[PlasmaDomain::be_x]+m_pd.m_eqs->grid(IdealMHD::bi_x),
                                        m_pd.m_grids[PlasmaDomain::be_y]+m_pd.m_eqs->grid(IdealMHD::bi_y))).abs();
    joule_heating = electrical_resistivity*current_density*current_density;
    grids_dt[4] += joule_heating;
}

void AnomalousResistivity::computeDiffusion(){
    const Grid &dx = m_pd.m_grids[PlasmaDomain::d_x], &dy = m_pd.m_grids[PlasmaDomain::d_y];
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);
    double rk_timestep = m_pd.epsilon*m_pd.m_eqs->getDT().min(m_pd.m_xl_dt,m_pd.m_yl_dt,m_pd.m_xu_dt,m_pd.m_yu_dt);
    curr_num_subcycles = (int)(1.0 + rk_timestep/(safety_factor*time_scale));
    diffusivity = one/(one/dx.square()+one/dy.square())/2./time_scale;
}

void AnomalousResistivity::computeTemplate(const Grid &b_x, const Grid &b_y){
    Grid b_mag_2d = (b_x.square() + b_y.square()).sqrt();
    Grid grad_b_frob_norm = ((m_pd.derivative1D(b_x,0)).square() + (m_pd.derivative1D(b_x,1)).square()
                        + (m_pd.derivative1D(b_y,0)).square() + (m_pd.derivative1D(b_y,1)).square()).sqrt();
    anomalous_template = (m_pd.laplacian(grad_b_frob_norm/b_mag_2d)).square();
    anomalous_template = (metric_coeff*anomalous_template).min(1.0);
    if(metric_smoothing){
        Grid smoothed_template = Grid::Zero(anomalous_template.rows(),anomalous_template.cols());
        for(int i=0; i<anomalous_template.rows(); i++){
            for(int j=0; j<anomalous_template.cols(); j++){
                for(int k=-kernel_radius; k<=kernel_radius; k++){
                    for(int l=-kernel_radius; l<=kernel_radius; l++){
                        if(i+k < 0 || j+l < 0 || i+k >= anomalous_template.rows() || j+l >= anomalous_template.cols()) continue;
                        smoothed_template(i,j) += smoothing_kernel(k+kernel_radius,l+kernel_radius)*anomalous_template(i+k,j+l);
                    }
                }
            }
        }
        anomalous_template = smoothed_template;
        anomalous_template = (10.0*anomalous_template).min(1.0);
        anomalous_template = anomalous_template.pow(1.5);
    }

}

std::string AnomalousResistivity::commandLineMessage() const
{
    return "Anomalous Resistivity Active";
}

void AnomalousResistivity::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids){
    if (output_to_file) {
        var_names.push_back("anomalous_diffusivity");
        var_grids.push_back(anomalous_template*diffusivity);
        var_names.push_back("anomalous_template");
        var_grids.push_back(anomalous_template);
        var_names.push_back("joule_heating");
        var_grids.push_back(avg_heating);
    }
}
