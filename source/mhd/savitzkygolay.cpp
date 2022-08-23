#include "savitzkygolay.hpp"
#include <string>
#include <cmath>
#include <cassert>

// constructor
SavitzkyGolay::SavitzkyGolay(std::string option,const Grid& grid)
{
    // initialize requested kernal and polynomial sizes
    switch (opt2ind(option)){
        case static_cast<int>(k55_p33): initialize_k55_p33(); break;
        case static_cast<int>(k33_p11): initialize_k33_p11(); break;
        default: assert(false);
    }
    m_Bx = (m_Kx-1)/2;
    m_By = (m_Ky-1)/2;
    // for given SG option, determine how many boundary cells and 
    identify_interior_edges(grid);
}

// return the index for the option corresponding to name
int SavitzkyGolay::opt2ind(std::string name)
{
    auto it = std::find(opts.begin(),opts.end(),name);
    if (it == opts.end()){
        std::cerr << name << "is not a valid option." << std::endl;
        assert(false);
    }
    return std::distance(opts.begin(),it);
}

// initialize coefficients for 5x5 kernel and 3x3 polynomial
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
    initialize_coeffs(coeffs,symmetry);
}

// initialize coefficients for 3x3 kernel and 1x1 polynomial
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
    initialize_coeffs(coeffs,symmetry);
}

// constructs full coefficient matrix from compacted coeffs and symmetry information
void SavitzkyGolay::initialize_coeffs(std::vector<std::vector<double>>& coeffs,std::vector<std::string>& symmetry)
{
    assert(coeffs.size()==symmetry.size());
    assert(coeffs.size() == ((m_Px+1)*(m_Py+1)));
    for (auto elem : coeffs) assert(elem.size() == (m_Kx*m_Ky+1)/2);
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

void SavitzkyGolay::identify_interior_edges(const Grid& grid)
{
    assert(grid.rows()>=m_Kx && grid.cols()>=m_Ky);
    m_grid_size_x = grid.rows();
    m_grid_size_y = grid.cols();
    m_xl = m_Bx;
    m_xu = grid.rows()-1-m_Bx;
    m_yl = m_By;
    m_yu = grid.cols()-1-m_By;
}

void SavitzkyGolay::process_boundary_cells(int ind_x, int ind_y, int& u_x, int& u_y, int& n_x, int& n_y) const
{
    // handle x boundary cells
    if (ind_x < m_xl){
        u_x = ind_x - m_xl;
        n_x = m_xl;
    }
    else if (ind_x > m_xu){
        u_x = ind_x - m_xu;
        n_x = m_xu;
    }
    else{
        u_x = 0;
        n_x = ind_x;
    }
    // handle y boundary cells
    if (ind_y < m_yl){
        u_y = ind_y - m_yl;
        n_y = m_yl;
    }
    else if (ind_y > m_yu){
        u_y = ind_y - m_yu;
        n_y = m_yu;
    }
    else{
        u_y = 0;
        n_y = ind_y;
    }
}

double SavitzkyGolay::smooth_cell_interior(const Grid& grid, int ind_x, int ind_y) const
{
    assert(ind_x >= m_xl && ind_x <= m_xu);
    assert(ind_y >= m_yl && ind_y <= m_yu);
    double result{0.};
    for (int kx : m_kx){
        for (int ky : m_ky){
            result += m_coeffs[0][0](kx,ky)*grid(ind_x+kx-m_Bx,ind_y+ky-m_By);
        }
    }
    return result;
}

double SavitzkyGolay::smooth_cell_boundary(const Grid& grid, int ind_x, int ind_y) const
{
    double result{0.};
    // determine boundary status
    bool is_bndry_x = (ind_x<m_xl || ind_x>m_xu) ? true : false;
    bool is_bndry_y = (ind_y<m_yl || ind_y>m_yu) ? true : false;
    assert(is_bndry_x || is_bndry_y);
    // get boundary cell positions relative to nearest interior cell
    int u_x, u_y; // cell positions relative to boundary
    int n_x, n_y; // nearest cell position to boundary cell
    process_boundary_cells(ind_x,ind_x,u_x,u_y,n_x,n_y);
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
            result += a_xy*x_contribution*y_contribution;
        }
    }
    return result;
}

double SavitzkyGolay::derivative_cell_interior(const Grid& grid, int ind_x, int ind_y, int mx, int my) const
{
    assert(ind_x >= m_xl && ind_x <= m_xu);
    assert(ind_y >= m_yl && ind_y <= m_yu);
    double result{0.};
    for (int kx : m_kx){
        for (int ky : m_ky){
            result += m_coeffs[mx][my](kx,ky)*grid(ind_x+kx-m_Bx,ind_y+ky-m_By);
        }
    }
    return result;
}

double SavitzkyGolay::derivative_cell_boundary(const Grid& grid, int ind_x, int ind_y, int mx, int my) const
{
    double result{0.};
    // determine boundary status
    bool is_bndry_x = (ind_x<m_xl || ind_x>m_xu) ? true : false;
    bool is_bndry_y = (ind_y<m_yl || ind_y>m_yu) ? true : false;
    assert(is_bndry_x || is_bndry_y);
    // get boundary cell positions relative to nearest interior cell
    int u_x, u_y; // cell positions relative to boundary
    int n_x, n_y; // nearest cell position to boundary cell
    process_boundary_cells(ind_x,ind_x,u_x,u_y,n_x,n_y);
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
            result += a_xy*x_contribution*y_contribution;
        }
    }
    return result;
}

// return smoothed grid - treats both interior and boundary cells
Grid SavitzkyGolay::smooth(const Grid& grid) const
{
    assert(grid.rows() == m_grid_size_x);
    assert(grid.cols() == m_grid_size_y);
    Grid interior = smooth_interior(grid);
    Grid boundary = smooth_boundary(grid);
    return interior + boundary;
}

Grid SavitzkyGolay::smooth(const Grid& grid, int xl, int xu, int yl, int yu) const
{
    assert(grid.rows() == m_grid_size_x);
    assert(grid.cols() == m_grid_size_y);
    assert(xl >= 0 && yl >= 0 && xu < grid.rows() && yu < grid.cols());
    Grid result = grid;
    #pragma omp parallel for
    for (int i=xl; i<=xu; i++){
        for (int j=yl; j<=yu; j++){
            if (i < m_xl || i > m_xu || j < m_yl || j > m_yu)
                result(i,j) = smooth_cell_boundary(grid,i,j);
            else 
                result(i,j) = smooth_cell_interior(grid,i,j);
        }
    }
    return result;
}

Grid SavitzkyGolay::smooth_interior(const Grid& grid) const
{
    assert(grid.rows() == m_grid_size_x);
    assert(grid.cols() == m_grid_size_y);
    // loop over interior cells
    Grid result = Grid::Zero(grid.rows(),grid.cols());
    #pragma omp parallel for
    for (int i=m_xl; i<=m_xu; i++){
        for (int j=m_yl; j<=m_yu; j++){
            result(i,j) = smooth_cell_interior(grid,i,j);
        }
    }
    return result;
}

Grid SavitzkyGolay::smooth_boundary(const Grid& grid) const
{
    assert(grid.rows() == m_grid_size_x);
    assert(grid.cols() == m_grid_size_y);
    Grid result {Grid::Zero(grid.rows(),grid.cols())};
    // lower x boundary
    for (int i=0; i<m_xl; i++){
        #pragma omp parallel for
        for (int j=0; j<grid.cols(); j++){
            result(i,j) = smooth_cell_boundary(grid,i,j);
        }
    }
    // upper x boundary
    for (int i=m_xu+1; i<grid.rows(); i++){
        #pragma omp parallel for
        for (int j=0; j<grid.cols(); j++){
            result(i,j) = smooth_cell_boundary(grid,i,j);
        }
    }
    // lower y boundary
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=0; j<m_yl; j++){
            result(i,j) = smooth_cell_boundary(grid,i,j);
        }
    }
    // upper y boundary
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=m_yu+1; j<grid.cols(); j++){
            result(i,j) = smooth_cell_boundary(grid,i,j);
        }
    }
    return result;
}

Grid SavitzkyGolay::derivative1D(const Grid& grid,int dim,double dr) const
{
    assert(grid.rows() == m_grid_size_x);
    assert(grid.cols() == m_grid_size_y);
    Grid interior = derivative1D_interior(grid,dim,dr);
    Grid boundary = derivative1D_boundary(grid,dim,dr);
    return interior + boundary;
}

Grid SavitzkyGolay::derivative1D_interior(const Grid& grid,int dim,double dr) const
{
    assert(grid.rows() == m_grid_size_x);
    assert(grid.cols() == m_grid_size_y);
    int mx = (dim==0) ? 1 : 0; // order of derivative along x
    int my = (dim==1) ? 1 : 0; // order of derivative along y
    // loop over interior cells
    Grid result = Grid::Zero(grid.rows(),grid.cols());
    #pragma omp parallel for
    for (int i=m_xl; i<=m_xu; i++){
        for (int j=m_yl; j<=m_yu; j++){
            result(i,j) = derivative_cell_interior(grid,i,j,mx,my);
        }
    }
    return result/dr;
}

Grid SavitzkyGolay::derivative1D_boundary(const Grid& grid,int dim,double dr) const
{
    assert(grid.rows() == m_grid_size_x);
    assert(grid.cols() == m_grid_size_y);
    int mx = (dim==0) ? 1 : 0; // order of derivative along x
    int my = (dim==1) ? 1 : 0; // order of derivative along y
    Grid result {Grid::Zero(grid.rows(),grid.cols())};
    // lower x boundary
    for (int i=0; i<m_xl; i++){
        #pragma omp parallel for
        for (int j=0; j<grid.cols(); j++){
            result(i,j) = derivative_cell_boundary(grid,i,j,mx,my);
        }
    }
    // upper x boundary
    for (int i=m_xu+1; i<grid.rows(); i++){
        #pragma omp parallel for
        for (int j=0; j<grid.cols(); j++){
            result(i,j) = derivative_cell_boundary(grid,i,j,mx,my);
        }
    }
    // lower y boundary
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=0; j<m_yl; j++){
            result(i,j) = derivative_cell_boundary(grid,i,j,mx,my);
        }
    }
    // upper y boundary
    #pragma omp parallel for
    for (int i=0; i<grid.rows(); i++){
        for (int j=m_yu+1; j<grid.cols(); j++){
            result(i,j) = derivative_cell_boundary(grid,i,j,mx,my);
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