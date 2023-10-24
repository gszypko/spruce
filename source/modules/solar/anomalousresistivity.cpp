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
    if(time_integrator == "") time_integrator = "euler";
    assert((time_integrator == "euler" || time_integrator == "rk2" || time_integrator == "rk4")
            && "Invalid time integrator given for Anomalous Resistivity module");
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
        else if(this_lhs == "time_integrator") time_integrator = this_rhs;
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

std::vector<Grid> AnomalousResistivity::computeTimeDerivatives(const Grid &bi_x, const Grid &bi_y, const Grid &bi_z, const Grid &be_x_lap, const Grid &be_y_lap, const Grid &be_z_lap){
    Grid b_x = m_pd.m_grids[PlasmaDomain::be_x]+bi_x;
    Grid b_y = m_pd.m_grids[PlasmaDomain::be_y]+bi_y;
    computeTemplate(b_x,b_y);
    Grid coeff = m_pd.m_ghost_zone_mask*anomalous_template*diffusivity;
    Grid electrical_resistivity = 4.0*PI/C/C*coeff;
    Grid current_density = C/(4.0*PI)*(m_pd.curl2D(b_x,b_y)).abs();
    Grid joule_heating = electrical_resistivity*current_density*current_density;
    return {coeff*(be_x_lap+m_pd.laplacian(bi_x)),
            coeff*(be_y_lap+m_pd.laplacian(bi_y)),
            coeff*(be_z_lap+m_pd.laplacian(bi_z)),
            joule_heating};
}

void AnomalousResistivity::iterateModule(double dt){
    computeDiffusion();
    Grid bi_x = m_pd.m_eqs->grid(IdealMHD::bi_x), bi_y = m_pd.m_eqs->grid(IdealMHD::bi_y), bi_z = m_pd.m_eqs->grid(IdealMHD::bi_z);
    Grid thermal_energy = m_pd.m_eqs->grid(IdealMHD::thermal_energy);
    if(output_to_file) avg_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    double dt_subcycle = (dt/(double)curr_num_subcycles);
    Grid old_thermal_energy = thermal_energy;
    Grid be_x_lap = m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_x]);
    Grid be_y_lap = m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_y]);
    Grid be_z_lap = m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_z]);
    if(time_integrator == "euler") {
        for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
            std::vector<Grid> derivatives = computeTimeDerivatives(bi_x,bi_y,bi_z,be_x_lap,be_y_lap,be_z_lap);
            bi_x += dt_subcycle*derivatives[0];
            bi_y += dt_subcycle*derivatives[1];
            bi_z += dt_subcycle*derivatives[2];
            thermal_energy += dt_subcycle*derivatives[3];
        }
    }
    else if(time_integrator == "rk2") {
        std::vector<Grid> bi_half_step, derivatives;
        for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
            derivatives = computeTimeDerivatives(bi_x,bi_y,bi_z,be_x_lap,be_y_lap,be_z_lap);
            bi_half_step = {bi_x,bi_y,bi_z};
            bi_half_step[0] += (0.5*dt_subcycle)*derivatives[0];
            bi_half_step[1] += (0.5*dt_subcycle)*derivatives[1];
            bi_half_step[2] += (0.5*dt_subcycle)*derivatives[2];

            derivatives = computeTimeDerivatives(bi_half_step[0],bi_half_step[1],bi_half_step[2],be_x_lap,be_y_lap,be_z_lap);
            bi_x += dt_subcycle*derivatives[0];
            bi_y += dt_subcycle*derivatives[1];
            bi_z += dt_subcycle*derivatives[2];
            thermal_energy += dt_subcycle*derivatives[3];
        }
    }
    else if(time_integrator == "rk4") {
        std::vector<Grid> intermediate_bi, k1, k2, k3, k4;
        for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
            k1 = computeTimeDerivatives(bi_x,bi_y,bi_z,be_x_lap,be_y_lap,be_z_lap);

            intermediate_bi = {bi_x + 0.5*dt_subcycle*k1[0],
                                bi_y + 0.5*dt_subcycle*k1[1],
                                bi_z + 0.5*dt_subcycle*k1[2]};
            k2 = computeTimeDerivatives(intermediate_bi[0],intermediate_bi[1],intermediate_bi[2],be_x_lap,be_y_lap,be_z_lap);

            intermediate_bi = {bi_x + 0.5*dt_subcycle*k2[0],
                                bi_y + 0.5*dt_subcycle*k2[1],
                                bi_z + 0.5*dt_subcycle*k2[2]};
            k3 = computeTimeDerivatives(intermediate_bi[0],intermediate_bi[1],intermediate_bi[2],be_x_lap,be_y_lap,be_z_lap);

            intermediate_bi = {bi_x + dt_subcycle*k3[0],
                                bi_y + dt_subcycle*k3[1],
                                bi_z + dt_subcycle*k3[2]};
            k4 = computeTimeDerivatives(intermediate_bi[0],intermediate_bi[1],intermediate_bi[2],be_x_lap,be_y_lap,be_z_lap);

            bi_x += dt_subcycle*(k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0])/6.0;
            bi_y += dt_subcycle*(k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1])/6.0;
            bi_z += dt_subcycle*(k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2])/6.0;
            thermal_energy += dt_subcycle*(k1[3] + 2.0*k2[3] + 2.0*k3[3] + k4[3])/6.0;
        }
    }
    if(output_to_file) avg_heating = (thermal_energy - old_thermal_energy)/dt;
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) = thermal_energy;
    m_pd.m_eqs->grid(IdealMHD::bi_x) = bi_x;
    m_pd.m_eqs->grid(IdealMHD::bi_y) = bi_y;
    m_pd.m_eqs->grid(IdealMHD::bi_z) = bi_z;
    m_pd.m_eqs->propagateChanges();
}


// //std::vector<Grid> grids_dt {d_rho_dt,d_mom_x_dt,d_mom_y_dt,d_mom_z_dt,d_thermal_energy_dt,d_bi_x_dt,d_bi_y_dt,d_bi_z_dt};
// void AnomalousResistivity::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt){
//     computeDiffusion();
//     computeTemplate(m_pd.m_grids[PlasmaDomain::be_x]+grids[IdealMHD::bi_x],m_pd.m_grids[PlasmaDomain::be_y]+grids[IdealMHD::bi_y]);
//     Grid coeff = m_pd.m_ghost_zone_mask*anomalous_template*diffusivity;
//     grids_dt[5] += coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_x])+m_pd.laplacian(grids[IdealMHD::bi_x]));
//     grids_dt[6] += coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_y])+m_pd.laplacian(grids[IdealMHD::bi_y]));
//     grids_dt[7] += coeff*(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_z])+m_pd.laplacian(grids[IdealMHD::bi_z]));
//     Grid electrical_resistivity = 4.0*PI/C/C*coeff;
//     Grid current_density = C/(4.0*PI)*(m_pd.curl2D(m_pd.m_grids[PlasmaDomain::be_x]+m_pd.m_eqs->grid(IdealMHD::bi_x),
//                                         m_pd.m_grids[PlasmaDomain::be_y]+m_pd.m_eqs->grid(IdealMHD::bi_y))).abs();
//     joule_heating = electrical_resistivity*current_density*current_density;
//     grids_dt[4] += joule_heating;
// }

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
    return "Anomalous Resistivity Subcycles: " + std::to_string(curr_num_subcycles);
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
