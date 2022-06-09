#include "sgfilter.hpp"
#include "plasmadomain.hpp"
#include "constants.hpp"
#include "idealmhd.hpp"
#include <iostream>

SGFilter::SGFilter(PlasmaDomain &pd): Module(pd) {
    assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
        "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
}

void SGFilter::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
    for(int i=0; i<lhs.size(); i++){
        if(lhs[i] == "filter_interval") filter_interval = std::stod(rhs[i]);
        else std::cerr << lhs[i] << " config not recognized.\n";
    }
}

void SGFilter::postIterateModule(double dt){
    if(filter_interval > 0 && m_pd.m_iter != 0 && m_pd.m_iter%filter_interval == 0){
        applyFilter();
    }
}

std::string SGFilter::commandLineMessage() const
{
    return "SG Filtering On";
}

// Applies two-dimensional Sovitzky-Golay filter to given Grid vector
// with 5x5 window and 3x3-order fitting polynomials
void SGFilter::applyFilter()
{
  #pragma omp parallel
  {
    #pragma omp for
    for(int varname : {IdealMHD::rho, IdealMHD::thermal_energy}){
      singleVarSavitzkyGolay(m_pd.grid(varname));
    }
  }
  m_pd.m_eqs->propagateChanges();
}

void SGFilter::singleVarSavitzkyGolay(Grid &grid)
{
  assert(N_GHOST == 2 && "S-G filtering implementation assumes two ghost zones");
  static const double coeff[] =
    {+7.346939E-03,	-2.938776E-02,	-4.163265E-02,	-2.938776E-02,	+7.346939E-03,
    -2.938776E-02,	+1.175510E-01,	+1.665306E-01,	+1.175510E-01,	-2.938776E-02,
    -4.163265E-02,	+1.665306E-01,	+2.359184E-01,  +1.665306E-01,  -4.163265E-02,
    -2.938776E-02,  +1.175510E-01,  +1.665306E-01,  +1.175510E-01,  -2.938776E-02,
    +7.346939E-03,  -2.938776E-02,  -4.163265E-02,  -2.938776E-02,  +7.346939E-03}; //from Chandra Sekhar, 2015
  int xdim = grid.rows();
  int ydim = grid.cols();
  Grid filtered = grid;
  for (int center_i = m_pd.m_xl; center_i <= m_pd.m_xu; center_i++){
    for(int center_j = m_pd.m_yl; center_j <= m_pd.m_yu; center_j++){
      int i[] = {center_i-2, center_i-1, center_i, center_i+1, center_i+2};
      int j[] = {center_j-2, center_j-1, center_j, center_j+1, center_j+2};
      if(m_pd.x_bound_1 == PlasmaDomain::BoundaryCondition::Periodic && m_pd.x_bound_2 == PlasmaDomain::BoundaryCondition::Periodic){
        i[0] = (i[0]+xdim)%xdim;
        i[1] = (i[1]+xdim)%xdim;
        i[3] = (i[3]+xdim)%xdim;
        i[4] = (i[4]+xdim)%xdim;
      }
      if(m_pd.y_bound_1 == PlasmaDomain::BoundaryCondition::Periodic && m_pd.y_bound_2 == PlasmaDomain::BoundaryCondition::Periodic){
        j[0] = (j[0]+ydim)%ydim;
        j[1] = (j[1]+ydim)%ydim;
        j[3] = (j[3]+ydim)%ydim;
        j[4] = (j[4]+ydim)%ydim;
      }
      double filtered_val = 0.0;
      for(int u = 0; u < 5; u++){
        for(int v = 0; v < 5; v++){
          int curr_i = i[u], curr_j = j[v];
          filtered_val = filtered_val + coeff[u*5 + v]*grid(curr_j,curr_j);
        }
      }
      filtered(center_i,center_j) = filtered_val;
    }
  }
  grid = filtered;
}