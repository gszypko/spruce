#include "savitzkygolay.hpp"
#include <string>
#include <cmath>
#include <cassert>

// constructor
SavitzkyGolay::SavitzkyGolay()
{
    initialize_k33_p11();
    m_Bx = (m_Kx-1)/2;
    m_By = (m_Ky-1)/2;
    void identify_boundary_cells();
}

// process 5x5 kernel, 3x3 polynomial coefficients
void SavitzkyGolay::initialize_k55_p33()
{
    // define polynomial order and kernel sizes
    m_Kx = 5; m_kx.resize(m_Kx);
    m_Ky = 5; m_ky.resize(m_Ky);
    m_Px = 3; m_px.resize(m_Px+1);
    m_Py = 3; m_py.resize(m_Py+1);
    for (int i=0; i<m_Kx; i++) m_kx[i] = i;
    for (int i=0; i<m_Ky; i++) m_ky[i] = i;
    for (int i=0; i<m_Px+1; i++) m_px[i] = i;
    for (int i=0; i<m_Py+1; i++) m_py[i] = i;
    // define coefficients - symmetry/antisymmetry still needs to be accounted for
    std::vector<std::vector<double>> coeffs;
    coeffs.push_back({+7.346939E-03,	-2.938776E-02,	-4.163265E-02,	-2.938776E-02,	+7.346939E-03,	-2.938776E-02,	+1.175510E-01,	+1.665306E-01,	+1.175510E-01,	-2.938776E-02,	-4.163265E-02,	+1.665306E-01,	+2.359184E-01}); // 0 0
	coeffs.push_back({-7.142857E-03,	+5.714286E-02,	+0.000000E+00,	-5.714286E-02,	+7.142857E-03,	+2.857143E-02,	-2.285714E-01,	+0.000000E+00,	+2.285714E-01,	-2.857143E-02,	+4.047619E-02,	-3.238095E-01,	+0.000000E+00}); // 0 1
	coeffs.push_back({-1.224490E-02,	+6.122449E-03,	+1.224490E-02,	+6.122449E-03,	-1.224490E-02,	+4.897959E-02,	-2.448980E-02,	-4.897959E-02,	-2.448980E-02,	+4.897959E-02,	+6.938776E-02,	-3.469388E-02,	-6.938776E-02}); // 0 2
	coeffs.push_back({+7.142857E-03,	-1.428571E-02,	+0.000000E+00,	+1.428571E-02,	-7.142857E-03,	-2.857143E-02,	+5.714286E-02,	+0.000000E+00,	-5.714286E-02,	+2.857143E-02,	-4.047619E-02,	+8.095238E-02,	+0.000000E+00}); // 0 3
	coeffs.push_back({-7.142857E-03,	+2.857143E-02,	+4.047619E-02,	+2.857143E-02,	-7.142857E-03,	+5.714286E-02,	-2.285714E-01,	-3.238095E-01,	-2.285714E-01,	+5.714286E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 1 0
	coeffs.push_back({+6.944444E-03,	-5.555556E-02,	+0.000000E+00,	+5.555556E-02,	-6.944444E-03,	-5.555556E-02,	+4.444444E-01,	+0.000000E+00,	-4.444444E-01,	+5.555556E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 1 1
	coeffs.push_back({+1.190476E-02,	-5.952381E-03,	-1.190476E-02,	-5.952381E-03,	+1.190476E-02,	-9.523810E-02,	+4.761905E-02,	+9.523810E-02,	+4.761905E-02,	-9.523810E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 1 2
	coeffs.push_back({-6.944444E-03,	+1.388889E-02,	+0.000000E+00,	-1.388889E-02,	+6.944444E-03,	+5.555556E-02,	-1.111111E-01,	+0.000000E+00,	+1.111111E-01,	-5.555556E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 1 3
	coeffs.push_back({-1.224490E-02,	+4.897959E-02,	+6.938776E-02,	+4.897959E-02,	-1.224490E-02,	+6.122449E-03,	-2.448980E-02,	-3.469388E-02,	-2.448980E-02,	+6.122449E-03,	+1.224490E-02,	-4.897959E-02,	-6.938776E-02}); // 2 0
	coeffs.push_back({+1.190476E-02,	-9.523810E-02,	+0.000000E+00,	+9.523810E-02,	-1.190476E-02,	-5.952381E-03,	+4.761905E-02,	+0.000000E+00,	-4.761905E-02,	+5.952381E-03,	-1.190476E-02,	+9.523810E-02,	+0.000000E+00}); // 2 1
	coeffs.push_back({+2.040816E-02,	-1.020408E-02,	-2.040816E-02,	-1.020408E-02,	+2.040816E-02,	-1.020408E-02,	+5.102041E-03,	+1.020408E-02,	+5.102041E-03,	-1.020408E-02,	-2.040816E-02,	+1.020408E-02,	+2.040816E-02}); // 2 2
	coeffs.push_back({-1.190476E-02,	+2.380952E-02,	+0.000000E+00,	-2.380952E-02,	+1.190476E-02,	+5.952381E-03,	-1.190476E-02,	+0.000000E+00,	+1.190476E-02,	-5.952381E-03,	+1.190476E-02,	-2.380952E-02,	+0.000000E+00}); // 2 3
	coeffs.push_back({+7.142857E-03,	-2.857143E-02,	-4.047619E-02,	-2.857143E-02,	+7.142857E-03,	-1.428571E-02,	+5.714286E-02,	+8.095238E-02,	+5.714286E-02,	-1.428571E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 3 0
	coeffs.push_back({-6.944444E-03,	+5.555556E-02,	+0.000000E+00,	-5.555556E-02,	+6.944444E-03,	+1.388889E-02,	-1.111111E-01,	+0.000000E+00,	+1.111111E-01,	-1.388889E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 3 1
	coeffs.push_back({-1.190476E-02,	+5.952381E-03,	+1.190476E-02,	+5.952381E-03,	-1.190476E-02,	+2.380952E-02,	-1.190476E-02,	-2.380952E-02,	-1.190476E-02,	+2.380952E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 3 2
	coeffs.push_back({+6.944444E-03,	-1.388889E-02,	+0.000000E+00,	+1.388889E-02,	-6.944444E-03,	-1.388889E-02,	+2.777778E-02,	+0.000000E+00,	-2.777778E-02,	+1.388889E-02,	+0.000000E+00,	+0.000000E+00,	+0.000000E+00}); // 3 3
    std::vector<std::string> symmetry {"S","A","S","A","A","S","A","S","S","A","S","A","A","S","A","S"};
    // initialize size of coefficient container
    m_coeffs.resize(m_Px+1);
    for (auto& row : m_coeffs) row.resize(m_Py+1,Grid::Zero(m_Kx,m_Ky));
    // process symmetry of coefficients
    for (int px : m_px){
        for (int py : m_py){
            int poly_ind = px*(m_Px+1)+py;
            for (int kx : m_kx){
                for (int ky : m_ky){
                    int kernel_ind = kx*m_Kx+ky;
                    double val;
                    if (kernel_ind < coeffs[poly_ind].size()) val = coeffs[poly_ind][kernel_ind];
                    else{
                        val = coeffs[poly_ind][m_Kx*m_Ky-1-kernel_ind];
                        if (symmetry[poly_ind]=="A") val *= (-1.);
                    }
                    m_coeffs[px][py](kx,ky) = val;
                }
            }
        }
        
    }
}

// process 5x5 kernel, 3x3 polynomial coefficients
void SavitzkyGolay::initialize_k33_p11()
{
    // define polynomial order and kernel sizes
    m_Kx = 3; m_kx.resize(m_Kx);
    m_Ky = 3; m_ky.resize(m_Ky);
    m_Px = 1; m_px.resize(m_Px+1);
    m_Py = 1; m_py.resize(m_Py+1);
    for (int i=0; i<m_Kx; i++) m_kx[i] = i;
    for (int i=0; i<m_Ky; i++) m_ky[i] = i;
    for (int i=0; i<m_Px+1; i++) m_px[i] = i;
    for (int i=0; i<m_Py+1; i++) m_py[i] = i;
    // define coefficients - symmetry/antisymmetry still needs to be accounted for
    std::vector<std::vector<double>> coeffs;
    coeffs.push_back({+1.111111E-01,	+1.111111E-01,	+1.111111E-01,	+1.111111E-01,	+1.111111E-01}); // 0 0
	coeffs.push_back({-1.666667E-01,	+0.000000E+00,	+1.666667E-01,	-1.666667E-01,	+0.000000E+00}); // 0 1
	coeffs.push_back({-1.666667E-01,	-1.666667E-01,	-1.666667E-01,	+0.000000E+00,	+0.000000E+00}); // 1 0
	coeffs.push_back({+2.500000E-01,	+0.000000E+00,	-2.500000E-01,	+0.000000E+00,	+0.000000E+00}); // 1 1
    std::vector<std::string> symmetry {"S","A","A","S"};
    // initialize size of coefficient container
    m_coeffs.resize(m_Px+1);
    for (auto& row : m_coeffs) row.resize(m_Py+1,Grid::Zero(m_Kx,m_Ky));
    // process symmetry of coefficients
    for (int px : m_px){
        for (int py : m_py){
            int poly_ind = px*(m_Px+1)+py;
            for (int kx : m_kx){
                for (int ky : m_ky){
                    int kernel_ind = kx*m_Kx+ky;
                    double val;
                    if (kernel_ind < coeffs[poly_ind].size()) val = coeffs[poly_ind][kernel_ind];
                    else{
                        val = coeffs[poly_ind][m_Kx*m_Ky-1-kernel_ind];
                        if (symmetry[poly_ind]=="A") val *= (-1.);
                    }
                    m_coeffs[px][py](kx,ky) = val;
                }
            }
        }
        
    }
}

Grid SavitzkyGolay::smoothing(const Grid& grid) const
{
    // identify nearest interior cells to boundary
    int Bx_l = m_Bx;
    int Bx_r = grid.rows()-1-m_Bx;
    int By_l = m_By;
    int By_r = grid.rows()-1-m_By;
    // begin differentiation loop
    Grid result {Grid::Zero(grid.rows(),grid.cols())};
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=0; j<grid.cols(); j++){
            // determine boundary status
            bool is_bndry_x = (i<Bx_l || i>Bx_r) ? true : false;
            bool is_bndry_y = (j<By_l || j>By_r) ? true : false;
            // handle interior cells
            if (!is_bndry_x && !is_bndry_y){
                for (int kx : m_kx){
                    for (int ky : m_ky){
                        result(i,j) += m_coeffs[0][0](kx,ky)*grid(i+kx-m_Bx,j+ky-m_By);
                    }
                }
            }
            // handle near-boundary cells
            else{
                // get boundary cell positions relative to nearest interior cell
                int u_x, u_y; // cell positions relative to boundary
                int n_x, n_y; // nearest cell position to boundary cell
                // handle x boundary cells
                if (i < Bx_l){
                    u_x = i - Bx_l;
                    n_x = Bx_l;
                }
                else if (i > Bx_r){
                    u_x = i - Bx_r;
                    n_x = Bx_r;
                }
                else{
                    u_x = 0;
                    n_x = i;
                }
                // handle y boundary cells
                if (j < By_l){
                    u_y = j - By_l;
                    n_y = By_l;
                }
                else if (j > By_r){
                    u_y = j - By_r;
                    n_y = By_r;
                }
                else{
                    u_y = 0;
                    n_y = j;
                }
                // loop over coefficients
                int mx = 0;
                int my = 0;
                for (int px = mx; px<m_Px; px++){
                    for (int py = my; py<m_Py; py++){
                        double a_xy{0.};
                        for (int kx : m_kx){
                            for (int ky : m_ky){
                                a_xy += m_coeffs[px][py](kx,ky)*grid(n_x+kx-m_Bx,n_y+ky-m_By);
                            }
                        }
                        double x_contribution = fac(px)*pow(u_x,px-mx)/fac(px-mx);
                        double y_contribution = fac(py)*pow(u_y,py-my)/fac(py-my);
                        result(i,j) += a_xy*x_contribution*y_contribution;
                    }
                }
                        
            }
        }
    }
    return result;
}

Grid SavitzkyGolay::boundary_smoothing(const Grid& grid) const
{
    // identify nearest interior cells to boundary
    int Bx_l = m_Bx;
    int Bx_r = grid.rows()-1-m_Bx;
    int By_l = m_By;
    int By_r = grid.rows()-1-m_By;
    // begin differentiation loop
    Grid result {Grid::Zero(grid.rows(),grid.cols())};
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=0; j<grid.cols(); j++){
            // determine boundary status
            bool is_bndry_x = (i<Bx_l || i>Bx_r) ? true : false;
            bool is_bndry_y = (j<By_l || j>By_r) ? true : false;
            // handle interior cells
            if (!is_bndry_x && !is_bndry_y){
                result(i,j) = grid(i,j);
            }
            // handle near-boundary cells
            else{
                // get boundary cell positions relative to nearest interior cell
                int u_x, u_y; // cell positions relative to boundary
                int n_x, n_y; // nearest cell position to boundary cell
                // handle x boundary cells
                if (i < Bx_l){
                    u_x = i - Bx_l;
                    n_x = Bx_l;
                }
                else if (i > Bx_r){
                    u_x = i - Bx_r;
                    n_x = Bx_r;
                }
                else{
                    u_x = 0;
                    n_x = i;
                }
                // handle y boundary cells
                if (j < By_l){
                    u_y = j - By_l;
                    n_y = By_l;
                }
                else if (j > By_r){
                    u_y = j - By_r;
                    n_y = By_r;
                }
                else{
                    u_y = 0;
                    n_y = j;
                }
                // loop over coefficients
                int mx = 0;
                int my = 0;
                for (int px = mx; px<m_Px; px++){
                    for (int py = my; py<m_Py; py++){
                        double a_xy{0.};
                        for (int kx : m_kx){
                            for (int ky : m_ky){
                                a_xy += m_coeffs[px][py](kx,ky)*grid(n_x+kx-m_Bx,n_y+ky-m_By);
                            }
                        }
                        double x_contribution = fac(px)*pow(u_x,px-mx)/fac(px-mx);
                        double y_contribution = fac(py)*pow(u_y,py-my)/fac(py-my);
                        result(i,j) += a_xy*x_contribution*y_contribution;
                    }
                }
                        
            }
        }
    }
    return result;
}

Grid SavitzkyGolay::derivative1D(const Grid& grid,int dim,double dr) const
{
    // process dimension being integrated
    int mx = (dim==0) ? 1 : 0; // order of derivative along x
    int my = (dim==1) ? 1 : 0; // order of derivative along y
    // identify nearest interior cells to boundary
    int Bx_l = m_Bx;
    int Bx_r = grid.rows()-1-m_Bx;
    int By_l = m_By;
    int By_r = grid.rows()-1-m_By;
    // begin differentiation loop
    Grid result {Grid::Zero(grid.rows(),grid.cols())};
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=0; j<grid.cols(); j++){
            // determine boundary status
            bool is_bndry_x = (i<Bx_l || i>Bx_r) ? true : false;
            bool is_bndry_y = (j<By_l || j>By_r) ? true : false;
            // handle interior cells
            if (!is_bndry_x && !is_bndry_y){
                for (int kx : m_kx){
                    for (int ky : m_ky){
                        result(i,j) += m_coeffs[mx][my](kx,ky)*grid(i+kx-m_Bx,j+ky-m_By);
                    }
                }
            }
            // handle near-boundary cells
            else{
                // get boundary cell positions relative to nearest interior cell
                int u_x, u_y; // cell positions relative to boundary
                int n_x, n_y; // nearest cell position to boundary cell
                // handle x boundary cells
                if (i < Bx_l){
                    u_x = i - Bx_l;
                    n_x = Bx_l;
                }
                else if (i > Bx_r){
                    u_x = i - Bx_r;
                    n_x = Bx_r;
                }
                else{
                    u_x = 0;
                    n_x = i;
                }
                // handle y boundary cells
                if (j < By_l){
                    u_y = j - By_l;
                    n_y = By_l;
                }
                else if (j > By_r){
                    u_y = j - By_r;
                    n_y = By_r;
                }
                else{
                    u_y = 0;
                    n_y = j;
                }
                // loop over coefficients
                for (int px = mx; px<m_Px; px++){
                    for (int py = my; py<m_Py; py++){
                        double a_xy{0.};
                        for (int kx : m_kx){
                            for (int ky : m_ky){
                                a_xy += m_coeffs[px][py](kx,ky)*grid(n_x+kx-m_Bx,n_y+ky-m_By);
                            }
                        }
                        double x_contribution = fac(px)*pow(u_x,px-mx)/fac(px-mx);
                        double y_contribution = fac(py)*pow(u_y,py-my)/fac(py-my);
                        result(i,j) += a_xy*x_contribution*y_contribution;
                    }
                }
                        
            }
        }
    }
    return result/dr;
}

int SavitzkyGolay::fac(int val) const
{
    if (val == 0) return 1;
    else{
        int result = 1;
        while (val > 0){
            result *= val;
            val -= 1;
        }
        return result;
    }
}