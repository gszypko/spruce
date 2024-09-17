#include "module.hpp"
#include "plasmadomain.hpp"
#include "anomalousresistivity.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include "solarutils.hpp"
#include "utils.hpp"
#include <sstream>
#include <iostream>
#include <cmath>

AnomalousResistivity::AnomalousResistivity(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    time_scale = -1.0;
}

void AnomalousResistivity::setupModule(){
    assert((resistivity_model != "time_scale" || time_scale > 0.0) && "Anomalous Resistivity time_scale must be specified, and positive, if run in time_scale mode");
    diffusivity = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    joule_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    anomalous_template = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    avg_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    cumulative_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    if(metric_smoothing){
        assert(smoothing_sigma > 0.0 && "smoothing_sigma must be a positive number");
        kernel_radius = nearbyint(4.0*smoothing_sigma);
        int kernel_size = 2*kernel_radius + 1;
        smoothing_kernel = SolarUtils::GaussianGrid(kernel_size,kernel_size,0.0,1.0/(2.0*PI*smoothing_sigma*smoothing_sigma),
                                                    smoothing_sigma,smoothing_sigma,kernel_radius,kernel_radius);
    }
    else smoothing_kernel = Grid::Zero(1,1);
    computeTemplate(m_pd.m_grids[PlasmaDomain::be_x]+m_pd.m_eqs->grid(IdealMHD::bi_x), m_pd.m_grids[PlasmaDomain::be_y]+m_pd.m_eqs->grid(IdealMHD::bi_y));
    if(time_integrator == "") time_integrator = "euler";
    assert((time_integrator == "euler" || time_integrator == "rk2" || time_integrator == "rk4")
            && "Invalid time integrator given for Anomalous Resistivity module");
    assert((template_mode == "frobenius" || template_mode == "flood_fill") 
            && "Invalid template mode given for Anomalous Resistivity module");
    assert((resistivity_model == "time_scale" || resistivity_model == "syntelis_19"
                || resistivity_model == "gudiksen_11" || resistivity_model == "ys_94")
            && "Invalid resistivity model given for Anomalous Resistivity module");
}

void AnomalousResistivity::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "time_scale") time_scale = std::stod(rhs[i]);
        else if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "multispecies_output_mode") multispecies_output_mode = (this_rhs == "true");
        else if(this_lhs == "frobenius_metric_coeff") frobenius_metric_coeff = std::stod(this_rhs);
        else if(this_lhs == "smoothing_sigma") smoothing_sigma = std::stod(this_rhs);
        else if(this_lhs == "safety_factor") safety_factor = std::stod(this_rhs);
        else if(this_lhs == "metric_smoothing") metric_smoothing = (this_rhs == "true");
        else if(this_lhs == "time_integrator") time_integrator = this_rhs;
        else if(this_lhs == "template_mode") template_mode = this_rhs;
        else if(this_lhs == "resistivity_model") resistivity_model = this_rhs;
        else if(this_lhs == "resistivity_model_params"){
            std::vector<std::string> params_str = splitString(this_rhs,',');
            std::vector<double> b;
            std::for_each(params_str.begin(),params_str.end(), [&b](const std::string &el) { b.push_back(std::stod(el)); });
            resistivity_model_params = b;
        }
        else if(this_lhs == "flood_fill_threshold") flood_fill_threshold = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

std::vector<Grid> AnomalousResistivity::computeTimeDerivatives(const Grid &bi_x, const Grid &bi_y, const Grid &bi_z, const Grid &be_x_lap, const Grid &be_y_lap, const Grid &be_z_lap){
    Grid b_x = m_pd.m_grids[PlasmaDomain::be_x]+bi_x;
    Grid b_y = m_pd.m_grids[PlasmaDomain::be_y]+bi_y;
    computeTemplate(b_x,b_y);
    computeDiffusion(bi_x,bi_x);
    Grid coeff = m_pd.m_ghost_zone_mask*anomalous_template*diffusivity;
    Grid electrical_resistivity = 4.0*PI/C/C*coeff;
    Grid current_density = C/(4.0*PI)*(m_pd.curl2D(bi_x,bi_y)).abs();
    Grid joule_heating = electrical_resistivity*current_density*current_density;
    return {coeff*(be_x_lap+m_pd.laplacian(bi_x)),
            coeff*(be_y_lap+m_pd.laplacian(bi_y)),
            coeff*(be_z_lap+m_pd.laplacian(bi_z)),
            joule_heating};
}

void AnomalousResistivity::preIterateModule(double dt){
    //Here, reproducing the check made in evolution.cpp to determine if a file output has just occured
    int old_time_iter = (int)((m_pd.m_time - dt)/m_pd.m_time_output_interval);
    int new_time_iter = (int)(m_pd.m_time/m_pd.m_time_output_interval);
    bool store_cond_1 = m_pd.m_iter_output_interval > 0 && m_pd.m_iter%m_pd.m_iter_output_interval == 0;
    bool store_cond_2 = m_pd.m_time_output_interval > 0.0 && new_time_iter > old_time_iter;
    if (store_cond_1 || store_cond_2) cumulative_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim); //If so, since the cumulative joule heating has been output, we reset the count to zero
}

void AnomalousResistivity::iterateModule(double dt){
    Grid bi_x = m_pd.m_eqs->grid(IdealMHD::bi_x), bi_y = m_pd.m_eqs->grid(IdealMHD::bi_y), bi_z = m_pd.m_eqs->grid(IdealMHD::bi_z);
    computeTemplate(m_pd.m_grids[PlasmaDomain::be_x]+bi_x,m_pd.m_grids[PlasmaDomain::be_y]+bi_y);
    computeDiffusion(bi_x,bi_x);
    computeNumSubcycles();
    Grid thermal_energy = m_pd.m_eqs->grid(IdealMHD::thermal_energy);
    if(output_to_file) {
        avg_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
        cumulative_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    }
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
    if(output_to_file){
        cumulative_heating += (thermal_energy - old_thermal_energy);
        avg_heating = (thermal_energy - old_thermal_energy)/dt;
    }
    m_pd.m_eqs->grid(IdealMHD::thermal_energy) = thermal_energy;
    m_pd.m_eqs->grid(IdealMHD::bi_x) = bi_x;
    m_pd.m_eqs->grid(IdealMHD::bi_y) = bi_y;
    m_pd.m_eqs->grid(IdealMHD::bi_z) = bi_z;
    m_pd.m_eqs->propagateChanges();
}

void AnomalousResistivity::computeNumSubcycles(){
    const Grid &dx = m_pd.m_grids[PlasmaDomain::d_x], &dy = m_pd.m_grids[PlasmaDomain::d_y];
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);

    // catch any instance of everywhere-zero resistivity
    if(!((anomalous_template*diffusivity).max() > 0.0)){
        curr_num_subcycles = 1;
        return;
    }

    // determine time scale of diffusive evolution, if not using a time scale to define the diffusivity
    if(resistivity_model != "time_scale")
        time_scale = (one/(one/dx.square()+one/dy.square())/2./(anomalous_template*diffusivity)).min(m_pd.m_xl_dt,m_pd.m_yl_dt,m_pd.m_xu_dt,m_pd.m_yu_dt);

    double rk_timestep = m_pd.epsilon*m_pd.m_eqs->getDT().min(m_pd.m_xl_dt,m_pd.m_yl_dt,m_pd.m_xu_dt,m_pd.m_yu_dt);
    curr_num_subcycles = (int)(1.0 + rk_timestep/(safety_factor*time_scale));
}

void AnomalousResistivity::computeDiffusion(const Grid &bi_x, const Grid &bi_y){
    const Grid &dx = m_pd.m_grids[PlasmaDomain::d_x], &dy = m_pd.m_grids[PlasmaDomain::d_y];
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);

    if(resistivity_model == "time_scale"){
        diffusivity = one/(one/dx.square()+one/dy.square())/2./time_scale;
    }
    else if(resistivity_model == "syntelis_19"){
        assert(resistivity_model_params.size() == 3 && "Syntelis-19 resistivity model requires three parameters: eta_0, eta_1, and J_crit");
        double eta_0 = resistivity_model_params[0], eta_1 = resistivity_model_params[1], J_crit = resistivity_model_params[2];
        Grid current_density = C/(4.0*PI)*(m_pd.curl2D(bi_x,bi_y)).abs();
        Grid current_scaling = current_density/J_crit;
        for(int i=0; i<current_scaling.rows(); i++) for(int j=0; j<current_scaling.cols(); j++) if(current_scaling(i,j) < 1.0) current_scaling(i,j) = 0.0;
        diffusivity = eta_1*current_scaling + eta_0;
    }
    else if(resistivity_model == "ys_94"){
        assert(resistivity_model_params.size() == 3 && "YS-94 resistivity model requires three parameters: v_c, alpha, and eta_max");
        double v_c = resistivity_model_params[0], alpha = resistivity_model_params[1], eta_max = resistivity_model_params[2];
        Grid current_density = C/(4.0*PI)*(m_pd.curl2D(bi_x,bi_y)).abs();
        Grid electron_drift = current_density/m_pd.m_eqs->grid(IdealMHD::n)/E_CHARGE;
        Grid vel_ratio = electron_drift/v_c;
        diffusivity = (alpha*((vel_ratio - 1).square())).min(eta_max);
        for(int i=0; i<diffusivity.rows(); i++) for(int j=0; j<diffusivity.cols(); j++) if(vel_ratio(i,j) <= 1.0) diffusivity(i,j) = 0.0;
    }
}

void AnomalousResistivity::computeTemplate(const Grid &b_x, const Grid &b_y){
    Grid b_mag_2d = (b_x.square() + b_y.square()).sqrt();
    if(template_mode == "frobenius"){
        Grid grad_b_frob_norm = ((m_pd.derivative1D(b_x,0)).square() + (m_pd.derivative1D(b_x,1)).square()
                            + (m_pd.derivative1D(b_y,0)).square() + (m_pd.derivative1D(b_y,1)).square()).sqrt();
        anomalous_template = (m_pd.laplacian(grad_b_frob_norm/b_mag_2d)).square();
        anomalous_template = (frobenius_metric_coeff*anomalous_template).min(1.0);
    } else {
        assert(template_mode == "flood_fill");
        std::vector<int> null_location = b_mag_2d.argmin(m_pd.m_xl,m_pd.m_yl,m_pd.m_xu,m_pd.m_yu);
        anomalous_template = b_mag_2d.floodFill({null_location},flood_fill_threshold);
    }
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
    }
    if(template_mode == "frobenius"){
        anomalous_template = (10.0*anomalous_template).min(1.0);
        anomalous_template = anomalous_template.pow(1.5);
    }
}

std::string AnomalousResistivity::commandLineMessage() const
{
    return "Anomalous Resistivity Subcycles: " + std::to_string(curr_num_subcycles);
}

void AnomalousResistivity::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids){
    if (multispecies_output_mode) {
        var_names.push_back("total_direct_heating");
        var_grids.push_back(cumulative_heating);
    } else if (output_to_file) {
        var_names.push_back("anomalous_diffusivity");
        var_grids.push_back(anomalous_template*diffusivity);
        var_names.push_back("anomalous_template");
        var_grids.push_back(anomalous_template);
        var_names.push_back("joule_heating");
        var_grids.push_back(avg_heating);
    }
}
