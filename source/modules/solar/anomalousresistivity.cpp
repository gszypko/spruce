#include "module.hpp"
#include "plasmadomain.hpp"
#include "anomalousresistivity.hpp"
#include "idealmhd.hpp"
#include "constants.hpp"
#include <sstream>
#include <iostream>

AnomalousResistivity::AnomalousResistivity(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
    time_scale = -1.0;
}

void AnomalousResistivity::setupModule(){
    assert(time_scale > 0.0 && "Anomalous Resistivity time_scale must be specified, and positive");
    anomalous_template = Grid::Ones(m_pd.m_xdim,m_pd.m_ydim);
    diffusivity = (m_pd.m_eqs->grid(IdealMHD::b_mag)
                    /time_scale/(m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_x]+m_pd.m_eqs->grid(IdealMHD::bi_x)).square()
                                + m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_y]+m_pd.m_eqs->grid(IdealMHD::bi_y)).square()
                                + m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_z]+m_pd.m_eqs->grid(IdealMHD::bi_z)).square()).sqrt()
                                ).min(m_pd.m_xl,m_pd.m_yl,m_pd.m_xu,m_pd.m_yu);
}

void AnomalousResistivity::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        std::string this_lhs = lhs[i];
        std::string this_rhs = rhs[i];
        if(lhs[i] == "time_scale") time_scale = std::stod(rhs[i]);
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

//std::vector<Grid> grids_dt {d_rho_dt,d_mom_x_dt,d_mom_y_dt,d_mom_z_dt,d_thermal_energy_dt,d_bi_x_dt,d_bi_y_dt,d_bi_z_dt};
void AnomalousResistivity::computeTimeDerivativesModule(const std::vector<Grid> &grids,std::vector<Grid> &grids_dt){
    computeTemplate();
    Grid coeff = m_pd.m_ghost_zone_mask*anomalous_template*diffusivity;
    grids_dt[5] += coeff*m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_x]+grids[IdealMHD::bi_x]);
    grids_dt[6] += coeff*m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_y]+grids[IdealMHD::bi_y]);
    grids_dt[7] += coeff*m_pd.laplacian(m_pd.m_grids[PlasmaDomain::be_z]+grids[IdealMHD::bi_z]);
}

void AnomalousResistivity::computeTemplate(){
}

std::string AnomalousResistivity::commandLineMessage() const
{
    std::ostringstream oss;
    oss.precision(4);
    oss << diffusivity;
    return "Anomalous Resistivity Active (diffusivity " + oss.str() + ")";
}

void AnomalousResistivity::fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids)
{
}
