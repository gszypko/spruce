#ifndef SAVITZKYGOLAY_HPP
#define SAVITZKYGOLAY_HPP

#include "grid.hpp"

class SavitzkyGolay
{
public:
    SavitzkyGolay();
    void initialize_k55_p33();
    void initialize_k33_p11();
    Grid boundary_smoothing(const Grid& grid) const;
    Grid derivative1D(const Grid& grid,int dim,double dr) const;
    int fac(int val) const;
private:
    std::vector<std::vector<Grid>> m_coeffs; // convolution coefficients: m_coeffs[px][py](kx,ky)
    int m_Kx, m_Ky; // kernel sizes along the respective axes
    int m_Px, m_Py; // polynomial order for the respective axes
    int m_Bx, m_By; // number of boundary cells along each axis
    std::vector<int> m_px, m_py; // valid indices for element accessing the vector m_coeffs
    std::vector<int> m_kx, m_ky; // valid indices for element accessing the Grid m_coeffs[px][py]
};

#endif