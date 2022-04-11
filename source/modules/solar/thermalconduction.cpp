//thermalconduction.hpp
//Header for the Thermal Conduction Module,
//an implementation of the abstract Module class
//Applies field-aligned Spitzer-Harm thermal conductivity
//with free-streaming saturation

#include "module.hpp"
#include "plasmadomain.hpp"
#include "thermalconduction.hpp"

ThermalConduction::ThermalConduction(PlasmaDomain &pd): Module(pd) {}

void ThermalConduction::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(this_lhs == "flux_saturation")  flux_saturation = (this_rhs == "true");
        else if(this_lhs == "epsilon")  epsilon = std::stod(this_rhs);
        else if(this_lhs == "dt_subcycle_min")  dt_subcycle_min = std::stod(this_rhs);
        else std::cerr << this_lhs << " config not recognized for Thermal Conduction Module.\n";
    }
}

void ThermalConduction::preIterateModule(double dt){
    curr_num_subcycles = numberSubcycles(dt);
}

void ThermalConduction::iterateModule(double dt){
    //Subcycle to simulate field-aligned thermal conduction
    Grid &m_thermal_energy = m_pd.m_grids[PlasmaDomain::thermal_energy], &m_temp = m_pd.m_grids[PlasmaDomain::temp], &m_press = m_pd.m_grids[PlasmaDomain::press];
    // Grid nonthermal_energy = 0.5*(m_grids[mom_x]*m_grids[v_x] + m_grids[mom_y]*m_grids[v_y]) + m_grids[mag_press];
    // Grid energy_floor = nonthermal_energy + thermal_energy_min; //To ensure non-negative thermal pressure
    Grid thermal_energy_next;
    for(int subcycle = 0; subcycle < curr_num_subcycles; subcycle++){
        Grid con_flux_x(m_pd.m_xdim,m_pd.m_ydim,0.0);
        Grid con_flux_y(m_pd.m_xdim,m_pd.m_ydim,0.0);
        fieldAlignedConductiveFlux(con_flux_x, con_flux_y, m_temp, m_pd.m_grids[PlasmaDomain::rho], m_pd.m_grids[PlasmaDomain::b_hat_x], m_pd.m_grids[PlasmaDomain::b_hat_y], KAPPA_0);
        if(flux_saturation) saturateConductiveFlux(con_flux_x, con_flux_y, m_pd.m_grids[PlasmaDomain::rho], m_temp);
        thermal_energy_next = m_thermal_energy - m_pd.m_ghost_zone_mask*(dt/(double)curr_num_subcycles)*(m_pd.derivative1D(con_flux_x,0)+m_pd.derivative1D(con_flux_y,1));
        m_thermal_energy = thermal_energy_next.max(m_pd.thermal_energy_min);
        m_pd.propagateChanges();
    }
}

//Computes number of subcycles necessary for the current iteration of the module
int ThermalConduction::numberSubcycles(double dt){
    Grid dt_subcycle;
    if(!flux_saturation){
        dt_subcycle = K_B/KAPPA_0*(m_pd.m_grids[PlasmaDomain::rho]/m_pd.m_ion_mass)*m_pd.m_grids[PlasmaDomain::d_x]*m_pd.m_grids[PlasmaDomain::d_y]/m_pd.m_grids[PlasmaDomain::temp].pow(2.5);
    } else {
        Grid kappa_modified(m_pd.m_xdim,m_pd.m_ydim,0.0);
        Grid con_flux_x(m_pd.m_xdim,m_pd.m_ydim,0.0);
        Grid con_flux_y(m_pd.m_xdim,m_pd.m_ydim,0.0);
        Grid field_temp_gradient = m_pd.derivative1D(m_pd.m_grids[PlasmaDomain::temp],0)*m_pd.m_grids[PlasmaDomain::b_hat_x]
                                + m_pd.derivative1D(m_pd.m_grids[PlasmaDomain::temp],1)*m_pd.m_grids[PlasmaDomain::b_hat_y];
        fieldAlignedConductiveFlux(con_flux_x, con_flux_y, m_pd.m_grids[PlasmaDomain::temp], m_pd.m_grids[PlasmaDomain::rho],
                                    m_pd.m_grids[PlasmaDomain::b_hat_x], m_pd.m_grids[PlasmaDomain::b_hat_y], KAPPA_0);
        if(flux_saturation) saturateConductiveFlux(con_flux_x, con_flux_y, m_pd.m_grids[PlasmaDomain::rho], m_pd.m_grids[PlasmaDomain::temp]);
        Grid flux_mag = (con_flux_x.square() + con_flux_y.square()).sqrt();
        kappa_modified = (flux_mag/field_temp_gradient).abs();
        if(field_temp_gradient.abs().max() == 0.0) return 0;
        else {
            dt_subcycle = K_B/kappa_modified*(m_pd.m_grids[PlasmaDomain::rho]/m_pd.m_ion_mass)*m_pd.m_grids[PlasmaDomain::d_x]*m_pd.m_grids[PlasmaDomain::d_y];
        }
    }
    double min_dt_subcycle = std::max(epsilon*dt_subcycle.min(m_pd.m_xl,m_pd.m_yl,m_pd.m_xu,m_pd.m_yu),dt_subcycle_min);
    return (int)(dt/min_dt_subcycle) + 1;
}

//Computes 1D cell-centered conductive flux from temperature "temp"
//Flux computed in direction indicated by "index": 0 for x, 1 for y
//k0 is conductive coefficient
Grid ThermalConduction::oneDimConductiveFlux(const Grid &temp, const Grid &rho, double k0, int index){
    Grid kappa_max = m_pd.m_grids[PlasmaDomain::d_x]*m_pd.m_grids[PlasmaDomain::d_y]*K_B*(rho/m_pd.m_ion_mass)/dt_subcycle_min;
    return -(k0*temp.pow(5.0/2.0)).min(kappa_max)*m_pd.derivative1D(temp,index);
}

//Computes cell-centered, field-aligned conductive flux from temperature "temp"
//temp is temperature Grid
//b_hat_x, b_hat_y are the components of the *unit* vector b_hat
//k0 is conductive coefficient
//Output is written to flux_out_x and flux_out_y
void ThermalConduction::fieldAlignedConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &temp, const Grid &rho,
                                    const Grid &b_hat_x, const Grid &b_hat_y, double k0){
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
void ThermalConduction::saturateConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &rho, const Grid &temp){
    Grid sat_flux_mag = (1.0/6.0)*(3.0/2.0)*(rho/m_pd.m_ion_mass)*(K_B*temp).pow(1.5)/std::sqrt(M_ELECTRON);
    Grid flux_mag = (flux_out_x.square() + flux_out_y.square()).sqrt();
    Grid scale_factor = sat_flux_mag /((sat_flux_mag.square() + flux_mag.square()).sqrt());
    flux_out_x *= scale_factor;
    flux_out_y *= scale_factor;
}

std::string ThermalConduction::commandLineMessage() const
{
    return "Thermal Subcycles: " + std::to_string(curr_num_subcycles);
}