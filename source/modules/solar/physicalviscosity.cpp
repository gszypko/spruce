#include "module.hpp"
#include "plasmadomain.hpp"
#include "physicalviscosity.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <iostream>
#include <cassert>

PhysicalViscosity::PhysicalViscosity(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    avg_heating = Grid::Zero(1,1);
    coeff_grid = Grid::Zero(1,1);
}

void PhysicalViscosity::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "output_to_file") output_to_file = (this_rhs == "true");
        else if(this_lhs == "coeff") coeff = std::stod(this_rhs);
        else if(this_lhs == "ramp_length") ramp_length = std::stod(this_rhs);
        else if(this_lhs == "buffer_length") buffer_length = std::stod(this_rhs);
        else if(this_lhs == "epsilon") epsilon = std::stod(this_rhs);
        else if(this_lhs == "heating_on") heating_on = (this_rhs == "true");
        else if(this_lhs == "force_on") force_on = (this_rhs == "true");
        else if(this_lhs == "inactive_mode") inactive_mode = (this_rhs == "true");
        else if(this_lhs == "gradient_correction") gradient_correction = (rhs[i] == "true");
        else if(this_lhs == "time_integrator") time_integrator = this_rhs;
        else if(this_lhs == "ms_electron_heating_fraction") ms_electron_heating_fraction = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized.\n";
    }
}

void PhysicalViscosity::setupModule(){
    coeff_grid = constructCoefficientGrid(coeff,ramp_length,buffer_length);
    avg_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    avg_force = std::vector<Grid>(3,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    num_subcycles = 1;
    if(time_integrator == "") time_integrator = "euler";
    assert((time_integrator == "euler" || time_integrator == "rk2")
            && "Invalid time integrator given for Physical Viscosity module");
    assert(ms_electron_heating_fraction >= 0.0 && ms_electron_heating_fraction <= 1.0 && "Physical Viscosity MS electron heating fraction must be between 0 and 1");
}

Grid PhysicalViscosity::computeHeating(const std::vector<Grid>& v, const std::vector<Grid>& b_hat, const Grid& temperature){
    if(coeff == 0.0) return Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    else {
        const Grid &b_hat_x = b_hat[0], &b_hat_y = b_hat[1], &b_hat_z = b_hat[2];        
        const Grid &v_x = v[0], &v_y = v[1], &v_z = v[2];
        Grid del_x_v_x = m_pd.derivative1D(v_x,0), del_x_v_y = m_pd.derivative1D(v_y,0), del_x_v_z = m_pd.derivative1D(v_z,0),
             del_y_v_x = m_pd.derivative1D(v_x,1), del_y_v_y = m_pd.derivative1D(v_y,1), del_y_v_z = m_pd.derivative1D(v_z,1);
        Grid heating = 3.0*coeff_grid*temperature.pow(2.5)*((
            b_hat_x*(b_hat_x*del_x_v_x + b_hat_y*del_y_v_x)
            + b_hat_y*(b_hat_x*del_x_v_y + b_hat_y*del_y_v_y)
            + b_hat_z*(b_hat_x*del_x_v_z + b_hat_y*del_y_v_z)
            - (del_x_v_x + del_y_v_y)/3.0).square());
        // heating = scalar*(
        //     (1.0/3.0 - b_hat_x.square())*del_x_v_x + (1.0/3.0 - b_hat_y.square())*del_y_v_y
        //     -b_hat_x*b_hat_y*del_y_v_x - b_hat_y*b_hat_x*del_x_v_y - b_hat_z*b_hat_x*del_x_v_z - b_hat_z*b_hat_y*del_y_v_z
        // );
        return m_pd.m_ghost_zone_mask*heating;
    }
}

int PhysicalViscosity::computeViscousSubcycles(double dt){
    Grid v_mag = (m_pd.m_eqs->grid(IdealMHD::v_x).square()+m_pd.m_eqs->grid(IdealMHD::v_y).square()).sqrt();
    Grid directional_factor = 3.0*(m_pd.m_eqs->grid(IdealMHD::b_hat_x)*m_pd.m_eqs->grid(IdealMHD::v_x)/v_mag
                                    + m_pd.m_eqs->grid(IdealMHD::b_hat_y)*m_pd.m_eqs->grid(IdealMHD::v_y)/v_mag).square()
                                + 1.0; //3cos^2(theta) + 1
    for(int i=0; i<m_pd.m_xdim; i++){
        for(int j=0; j<m_pd.m_ydim; j++){
            if(v_mag(i,j) == 0.0) directional_factor(i,j) = 4.0; //set to maximum to be safe, in the ambiguous local v=0 case
        }
    }
    Grid timescale = (m_pd.m_grids[PlasmaDomain::d_x]*m_pd.m_grids[PlasmaDomain::d_y])*m_pd.m_eqs->grid(IdealMHD::rho)/(directional_factor*3.0*coeff_grid.max(0.001*coeff)*(m_pd.m_eqs->grid(IdealMHD::temp).pow(2.5)));
    double min_dt_subcycle = epsilon*timescale.min(m_pd.m_xl,m_pd.m_yl,m_pd.m_xu,m_pd.m_yu);
    return (int)(dt/min_dt_subcycle) + 1;
}

std::vector<Grid> PhysicalViscosity::computeViscousForce(const std::vector<Grid>& v, const std::vector<Grid>& b_hat, const Grid& temperature){
    std::vector<Grid> result(3,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    double delta;
    // std::vector<int> b_idxs = {IdealMHD::b_hat_x,IdealMHD::b_hat_y,IdealMHD::b_hat_z};
    const Grid &b_hat_x = b_hat[0], &b_hat_y = b_hat[1], &b_hat_z = b_hat[2];        
    const Grid &v_x = v[0], &v_y = v[1], &v_z = v[2];
    Grid del_x_v_x = m_pd.derivative1D(v_x,0), del_x_v_y = m_pd.derivative1D(v_y,0), del_x_v_z = m_pd.derivative1D(v_z,0),
            del_y_v_x = m_pd.derivative1D(v_x,1), del_y_v_y = m_pd.derivative1D(v_y,1), del_y_v_z = m_pd.derivative1D(v_z,1);
    Grid del_x_b_hat_y = m_pd.derivative1D(b_hat_y,0), del_y_b_hat_x = m_pd.derivative1D(b_hat_x,1),
            del_x_b_hat_x = m_pd.derivative1D(b_hat_x,0), del_y_b_hat_y = m_pd.derivative1D(b_hat_y,1);
    Grid same_factor_off_diag = (b_hat_x*(b_hat_y*del_y_v_x)
                                + b_hat_y*(b_hat_x*del_x_v_y)
                                + b_hat_z*(b_hat_x*del_x_v_z + b_hat_y*del_y_v_z)
                                - (del_x_v_x + del_y_v_y)/3.0);
    Grid same_factor_diag = (b_hat_x*(b_hat_x*del_x_v_x )
                                + b_hat_y*( b_hat_y*del_y_v_y));
    std::vector<Grid> grad_b_terms_diag = {
        2.0*b_hat_x*del_x_b_hat_x*del_x_v_x + b_hat_x*b_hat_x*m_pd.secondDerivative1D(v_x,0)
        + 2.0*b_hat_y*del_x_b_hat_y*del_y_v_y + b_hat_y*b_hat_y*m_pd.derivative1D(del_y_v_y,0),
        2.0*b_hat_x*del_y_b_hat_x*del_x_v_x + b_hat_x*b_hat_x*m_pd.derivative1D(del_x_v_x,1)
        + 2.0*b_hat_y*del_y_b_hat_y*del_y_v_y + b_hat_y*b_hat_y*m_pd.secondDerivative1D(v_y,1)
    }; // = divergence of same_factor_diag
    for(int j=0; j<3; j++){
        for(int i=0; i<2; i++){
            if(i == j) delta = 1.0/3.0;
            else delta = 0.0;
            //NOTE: Calculation is split up into three terms to expose the direct second-order velocity gradients (in grad_b_terms_diag)
            // result[j] += grad(b terms) dot (tensor)
            // result[j] += (b terms) div (tensor)
            if(gradient_correction){
                result[j] += m_pd.derivative1D(
                    temperature.pow(2.5)*coeff_grid
                    *(delta - b_hat[i]*b_hat[j])
                    *same_factor_off_diag, i);
                result[j] += grad_b_terms_diag[i]*temperature.pow(2.5)*coeff_grid
                    *(delta - b_hat[i]*b_hat[j]);
                result[j] += same_factor_diag * m_pd.derivative1D(
                    temperature.pow(2.5)*coeff_grid
                    *(delta - b_hat[i]*b_hat[j]), i);
            } else {
                result[j] += m_pd.derivative1D(
                    temperature.pow(2.5)
                    *(delta - b_hat[i]*b_hat[j])
                    *same_factor_off_diag, i);
                result[j] += grad_b_terms_diag[i]*temperature.pow(2.5)
                    *(delta - b_hat[i]*b_hat[j]);
                result[j] += same_factor_diag * m_pd.derivative1D(
                    temperature.pow(2.5)
                    *(delta - b_hat[i]*b_hat[j]), i);
            }
        }
        if(gradient_correction){
            result[j] *= -3.0*m_pd.m_ghost_zone_mask;
        } else{
            result[j] *= -3.0*coeff_grid*m_pd.m_ghost_zone_mask;
        }
    }
    return result;
}

void PhysicalViscosity::preIterateModule(double dt){
    if(heating_on) avg_heating = computeHeating({m_pd.m_eqs->grid(IdealMHD::v_x),m_pd.m_eqs->grid(IdealMHD::v_y),m_pd.m_eqs->grid(IdealMHD::v_z)},
                                        {m_pd.m_eqs->grid(IdealMHD::b_hat_x),m_pd.m_eqs->grid(IdealMHD::b_hat_y),m_pd.m_eqs->grid(IdealMHD::b_hat_z)},
                                        m_pd.m_eqs->grid(IdealMHD::temp));
}

void PhysicalViscosity::iterateModule(double dt){
    num_subcycles = computeViscousSubcycles(dt);
    double dt_subcycle = dt/num_subcycles;
    avg_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
    avg_force = std::vector<Grid>(3,Grid::Zero(m_pd.m_xdim,m_pd.m_ydim));
    Grid heating;
    std::vector<Grid> force;

    Grid v_x = m_pd.m_eqs->grid(IdealMHD::v_x), v_y = m_pd.m_eqs->grid(IdealMHD::v_y), v_z = m_pd.m_eqs->grid(IdealMHD::v_z);
    Grid mom_x = m_pd.m_eqs->grid(IdealMHD::mom_x), mom_y = m_pd.m_eqs->grid(IdealMHD::mom_y), mom_z = m_pd.m_eqs->grid(IdealMHD::mom_z);
    Grid temperature = m_pd.m_eqs->grid(IdealMHD::temp), thermal_energy = m_pd.m_eqs->grid(IdealMHD::thermal_energy);

    const Grid& b_hat_x = m_pd.m_eqs->grid(IdealMHD::b_hat_x), b_hat_y = m_pd.m_eqs->grid(IdealMHD::b_hat_y), b_hat_z = m_pd.m_eqs->grid(IdealMHD::b_hat_z);
    const Grid& m_n = m_pd.m_eqs->grid(IdealMHD::n), rho = m_pd.m_eqs->grid(IdealMHD::rho);

    if(time_integrator == "euler") for(int i=0; i<num_subcycles; i++){
        if(heating_on){
            heating = computeHeating({v_x,v_y,v_z},{b_hat_x,b_hat_y,b_hat_z},temperature);
            avg_heating += heating/num_subcycles;
        }
        if(force_on){
            force = computeViscousForce({v_x,v_y,v_z},{b_hat_x,b_hat_y,b_hat_z},temperature);
            for(int k=0; k<3; k++) avg_force[k] += force[k]/num_subcycles;
        }
        if(heating_on && !inactive_mode){
            thermal_energy += heating*dt_subcycle;
            if(m_pd.m_multispecies_mode){
                if(ms_electron_heating_fraction < 1.0) m_pd.m_cumulative_ion_heating += (1.0 - ms_electron_heating_fraction)*heating*dt_subcycle;
                if(ms_electron_heating_fraction > 0.0) m_pd.m_cumulative_electron_heating += ms_electron_heating_fraction*heating*dt_subcycle;
            }

            thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
            temperature = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*m_n)).max(m_pd.temp_min);
        }
        if(force_on && !inactive_mode){
            mom_x += force[0]*dt_subcycle;
            mom_y += force[1]*dt_subcycle;
            mom_z += force[2]*dt_subcycle;
            v_x = mom_x / rho;
            v_y = mom_y / rho;
            v_z = mom_z / rho;
        }
    }
    else if(time_integrator == "rk2"){
        Grid half_step_thermal = thermal_energy, half_step_temp = temperature;
        Grid v_x_half = v_x, v_y_half = v_y, v_z_half = v_z, mom_x_half = mom_x, mom_y_half = mom_y, mom_z_half = mom_z;
        for(int i=0; i<num_subcycles; i++){
            // start with half-step
            if(heating_on){
                heating = computeHeating({v_x,v_y,v_z},{b_hat_x,b_hat_y,b_hat_z},temperature);
            }
            if(force_on){
                force = computeViscousForce({v_x,v_y,v_z},{b_hat_x,b_hat_y,b_hat_z},temperature);
            }
            if(heating_on && !inactive_mode){
                half_step_thermal = thermal_energy + heating*0.5*dt_subcycle;
                half_step_thermal = half_step_thermal.max(m_pd.thermal_energy_min);
                half_step_temp = ((m_pd.m_adiabatic_index - 1.0)*half_step_thermal/(2*K_B*m_n)).max(m_pd.temp_min);
            }
            if(force_on && !inactive_mode){
                mom_x_half = mom_x + force[0]*0.5*dt_subcycle;
                mom_y_half = mom_y + force[1]*0.5*dt_subcycle;
                mom_z_half = mom_z + force[2]*0.5*dt_subcycle;
                v_x_half = mom_x_half / rho;
                v_y_half = mom_y_half / rho;
                v_z_half = mom_z_half / rho;
            }
            //now compute full time-step
            if(heating_on){
                heating = computeHeating({v_x_half,v_y_half,v_z_half},{b_hat_x,b_hat_y,b_hat_z},half_step_temp);
                avg_heating += heating/num_subcycles;
            }
            if(force_on){
                force = computeViscousForce({v_x_half,v_y_half,v_z_half},{b_hat_x,b_hat_y,b_hat_z},half_step_temp);
                for(int k=0; k<3; k++) avg_force[k] += force[k]/num_subcycles;
            }
            if(heating_on && !inactive_mode){
                thermal_energy += heating*dt_subcycle;
                if(m_pd.m_multispecies_mode) {
                    if(ms_electron_heating_fraction < 1.0) m_pd.m_cumulative_ion_heating += (1.0 - ms_electron_heating_fraction)*heating*dt_subcycle;
                    if(ms_electron_heating_fraction > 0.0) m_pd.m_cumulative_electron_heating += ms_electron_heating_fraction*heating*dt_subcycle;
                }

                thermal_energy = thermal_energy.max(m_pd.thermal_energy_min);
                temperature = ((m_pd.m_adiabatic_index - 1.0)*thermal_energy/(2*K_B*m_n)).max(m_pd.temp_min);
            }
            if(force_on && !inactive_mode){
                mom_x += force[0]*dt_subcycle;
                mom_y += force[1]*dt_subcycle;
                mom_z += force[2]*dt_subcycle;
                v_x = mom_x / rho;
                v_y = mom_y / rho;
                v_z = mom_z / rho;
            }
        }
    }

    m_pd.m_eqs->grid(IdealMHD::thermal_energy) = thermal_energy;
    m_pd.m_eqs->grid(IdealMHD::mom_x) = mom_x;
    m_pd.m_eqs->grid(IdealMHD::mom_y) = mom_y;
    m_pd.m_eqs->grid(IdealMHD::mom_z) = mom_z;
    m_pd.m_eqs->propagateChanges();    
}

Grid PhysicalViscosity::constructCoefficientGrid(double strength,double ramp_length,double buffer_length) const
{
    if(ramp_length == 0.0) return strength*Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);
    // initialize grids and references to PlasmaDomain grids
    const Grid& x = m_pd.grid("pos_x");
    const Grid& y = m_pd.grid("pos_y");
    // note: the choice of the 2.3 coefficient ensures that the exponential reaches about 1% after 1 full length
    // left boundary
    Grid result = 1.01*(-2.3*((-x + (x.min()+ramp_length+buffer_length)).max(0.0)/ramp_length).square()).exp() - 0.01;
    // right boundary
    result = result.min(1.01*(-2.3*((x - (x.max()-ramp_length-buffer_length)).max(0.0)/ramp_length).square()).exp() - 0.01);
    // top boundary
    result = result.min(1.01*(-2.3*((y - (y.max()-ramp_length-buffer_length)).max(0.0)/ramp_length).square()).exp() - 0.01);
    // bottom boundary
    result = result.min(1.01*(-2.3*((-y + (y.min()+ramp_length+buffer_length)).max(0.0)/ramp_length).square()).exp() - 0.01);
    // piecewise minimum calls ensure that the corners a constructed correctly
    result = result.max(0.0);
    return strength*result;
}


std::string PhysicalViscosity::commandLineMessage() const
{
    std::string message = "";
    if(heating_on){
        message += "Viscous Heating";
        if(coeff == 0.0) message += " Zero";
        else message += " On";
        if(force_on) message += ", ";
    }
    if(force_on){
        message += "Viscous Force";
        if(coeff == 0.0) message += " Zero";
        else message += " On";
    }
    if(force_on || heating_on) message += ", " + std::to_string(num_subcycles) + " Subcycle(s)";
    if(inactive_mode) message += " (Not Applied)";
    return message;
}

void PhysicalViscosity::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
    if(output_to_file){
        if(avg_heating.size() == 1) avg_heating = Grid::Zero(m_pd.m_xdim,m_pd.m_ydim);
        if(heating_on){ 
            var_names.push_back("viscous_heating");
            var_grids.push_back(avg_heating);
        }
        if(force_on){
            var_names.push_back("viscous_force_x");
            var_grids.push_back(avg_force[0]);
            var_names.push_back("viscous_force_y");
            var_grids.push_back(avg_force[1]);
            var_names.push_back("viscous_force_z");
            var_grids.push_back(avg_force[2]);
        }
    }
}
