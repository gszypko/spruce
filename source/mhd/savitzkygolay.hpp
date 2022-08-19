#ifndef SAVITZKYGOLAY_HPP
#define SAVITZKYGOLAY_HPP

#include "grid.hpp"
#include <vector>
#include <string>

// compacted coefficients from:
    // publication - http://aip.scitation.org/doi/abs/10.1063/1.4940262
    // repository - https://sites.google.com/site/chandraacads/resources/sg-filter/db

class SavitzkyGolay
{
public:
    // options for initialization
    static int opt2ind(std::string name);
    static inline std::vector<std::string> opts {"k55_p33","k33_p11"};
    enum opts {k55_p33,k33_p11};

    // construction
    SavitzkyGolay(std::string option = "k33_p11", const Grid& grid = Grid::Zero(4,4));

    // usage
    Grid smooth(const Grid& grid) const;
    Grid smooth_interior(const Grid& grid) const;
    Grid smooth_boundary(const Grid& grid) const;
    Grid derivative1D(const Grid& grid,int dim,double dr) const;
    Grid derivative1D_interior(const Grid& grid,int dim,double dr) const;
    Grid derivative1D_boundary(const Grid& grid,int dim,double dr) const;
private:
    // private members
    std::vector<std::vector<Grid>> m_coeffs; // convolution coefficients: m_coeffs[px][py](kx,ky)
    int m_Kx, m_Ky; // kernel sizes along the respective axes
    int m_Px, m_Py; // polynomial order for the respective axes
    int m_Bx, m_By; // number of boundary cells along each axis
    std::vector<int> m_px, m_py; // valid indices for element accessing the vector m_coeffs
    std::vector<int> m_kx, m_ky; // valid indices for element accessing the Grid m_coeffs[px][py]
    int m_grid_size_x, m_grid_size_y; // intended grid size for operations
    int m_xl, m_xu, m_yl, m_yu; // edges of interior domain

    // initialization
    void initialize_k55_p33();
    void initialize_k33_p11();
    void initialize_coeffs(std::vector<std::vector<double>>& coeffs,std::vector<std::string>& symmetry);

    // other functions
    int fac(int val) const;
    void identify_interior_edges(const Grid& grid);
    void process_boundary_cells(int ind_x, int ind_y, int& u_x, int& u_y, int& n_x, int& n_y) const;
    double smooth_cell_interior(const Grid& grid, int ind_x, int ind_y) const;
    double smooth_cell_boundary(const Grid& grid, int ind_x, int ind_y) const;
    double derivative_cell_interior(const Grid& grid, int ind_x, int ind_y, int mx, int my) const;
    double derivative_cell_boundary(const Grid& grid, int ind_x, int ind_y, int mx, int my) const;
};

#endif