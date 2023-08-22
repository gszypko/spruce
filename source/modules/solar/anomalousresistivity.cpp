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
        else if(this_lhs == "metric_smoothing") metric_smoothing = (this_rhs == "true");
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

//std::vector<Grid> grids_dt {d_rho_dt,d_mom_x_dt,d_mom_y_dt,d_mom_z_dt,d_thermal_energy_dt,d_bi_x_dt,d_bi_y_dt,d_bi_z_dt};
void AnomalousResistivity::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt){
    computeDiffusion(grids);
    computeTemplate(grids);
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

void AnomalousResistivity::computeDiffusion(const std::vector<Grid> &grids){
    const Grid &dx = m_pd.m_grids[PlasmaDomain::d_x], &dy = m_pd.m_grids[PlasmaDomain::d_y];
    Grid one = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);
    double rk_timestep = m_pd.epsilon*m_pd.m_eqs->getDT().min(m_pd.m_xl_dt,m_pd.m_yl_dt,m_pd.m_xu_dt,m_pd.m_yu_dt);
    if(rk_timestep > time_scale){
        std::cout << "Fluid time step exceeded anomalous diffusion time scale; using fluid time step as diffusion time scale\n";
        diffusivity = one/(one/dx.square()+one/dy.square())/2./(rk_timestep);
    } else {
        diffusivity = one/(one/dx.square()+one/dy.square())/2./time_scale;
    }
}

void AnomalousResistivity::computeTemplate(const std::vector<Grid> &grids){
    Grid b_mag_2d = ((m_pd.m_grids[PlasmaDomain::be_x] + grids[IdealMHD::bi_x]).square()
                    + (m_pd.m_grids[PlasmaDomain::be_y] + grids[IdealMHD::bi_y]).square()).sqrt();
    Grid grad_b_frob_norm = ((m_pd.derivative1D(m_pd.m_grids[PlasmaDomain::be_x] + grids[IdealMHD::bi_x],0)).square()
                        + (m_pd.derivative1D(m_pd.m_grids[PlasmaDomain::be_x] + grids[IdealMHD::bi_x],1)).square()
                        + (m_pd.derivative1D(m_pd.m_grids[PlasmaDomain::be_y] + grids[IdealMHD::bi_y],0)).square()
                        + (m_pd.derivative1D(m_pd.m_grids[PlasmaDomain::be_y] + grids[IdealMHD::bi_y],1)).square()).sqrt();
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
        var_grids.push_back(joule_heating);
    }
}
