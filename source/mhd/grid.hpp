//Class definition for 2D grid of double precision data
//Supports component-wise arithmetic, among other functions
//Built in usage of OpenMP parallelization
#ifndef GRID_HPP
#define GRID_HPP

#include <vector>
#include <string>
#include <functional>
#include <cassert>
#include <omp.h>
#include <sstream>
#include <iostream>

class Grid
{
public:
  // *** Construction
  Grid(size_t rows, size_t cols);
  Grid(size_t rows, size_t cols, double val);
  Grid(size_t rows, size_t cols, std::vector<double> data);
  Grid(std::vector<double> vec);
  Grid();
  static Grid Zero(size_t rows, size_t cols);
  static Grid Ones(size_t rows, size_t cols);
  static Grid Linspace(double start,double end,int num,int dim=1);
  // *** Element Accession
  double& operator()(size_t i, size_t j);
  double operator()(size_t i, size_t j) const;
  double& operator()(size_t el);
  double operator()(size_t el) const;
  // *** Basic Getters
  int rows() const;
  int cols() const;
  int size() const;
  std::vector<double> data() const;
  Grid col(size_t ind) const;
  Grid row(size_t ind) const;
  // *** Other Usage
  double min() const;
  double min(int il, int jl, int iu, int ju) const;
  Grid min(double b) const;
  Grid min(const Grid& b) const;
  double max() const;
  double max(int il, int jl, int iu, int ju) const;
  Grid max(double b) const;
  Grid max(const Grid& b) const;
  Grid square() const;
  Grid abs() const;
  Grid sqrt() const;
  Grid pow(double power) const;
  Grid log() const;
  Grid for_each(const Grid& grid,const std::function<double(double,double)>& fun) const;
  Grid for_each(double a,const std::function<double(double,double)>& fun) const;
  // *** New Math
  Grid vectorise(const Grid& grid,int dim=1) const;
  static Grid trapz(const Grid& x,const Grid& y);
  static Grid trapz2D(const Grid& x,const Grid& y,const Grid& z);
  static Grid trapzcum(const Grid& x, const Grid& y);
  static Grid bin_as_list(const Grid& indep,const Grid& dep, const Grid& bin_edges);
  static Grid bin_as_list(const Grid& indep,const Grid& dep, double N, Grid& bins);
  static Grid interp_as_list(const Grid& indep,const Grid& dep,const Grid& query);
  // *** Vector Products
  static Grid DotProduct2D(const std::vector<Grid>& a, const std::vector<Grid>& b);//Dot product of two 2D vectors
  static Grid CrossProduct2D(const std::vector<Grid>& a, const std::vector<Grid>& b);//Cross product of two vectors in xy-plane (result in z-direction)
  static std::vector<Grid> CrossProductZ2D(const Grid& a_z, const std::vector<Grid>& b);//Cross product of vector in z-direction with vector in xy-plane (result in xy-plane)
  static std::vector<Grid> CrossProduct2DZ(const std::vector<Grid>& a, const Grid& b_z);//Cross product of vector in xy-plane with vector in z-direction (result in xy-plane)
  static std::vector<Grid> CrossProduct(const std::vector<Grid>& a, const std::vector<Grid>& b); //General cross product between two 3D vector quantities
  // *** OPERATOR OVERRIDES
  Grid& operator+=(const Grid& b);
  Grid& operator+=(double b);
  Grid& operator-=(const Grid& b);
  Grid& operator-=(double b);
  Grid& operator*=(const Grid& b);
  Grid& operator*=(double b);
  Grid& operator/=(const Grid& b);
  Grid& operator/=(double b);
  friend std::ostream& operator<<(std::ostream& os, const Grid& a);
  friend Grid operator+(const Grid& a, const Grid& b);
  friend Grid operator+(const Grid& a, double b);
  friend Grid operator+(double a, const Grid& b);
  friend Grid operator-(const Grid& a, const Grid& b);
  friend Grid operator-(const Grid& a, double b);
  friend Grid operator-(double a, const Grid& b);
  friend Grid operator-(const Grid& a);
  friend Grid operator*(const Grid& a, const Grid& b);
  friend Grid operator*(const Grid& a, double b);
  friend Grid operator*(double a, const Grid& b);
  friend Grid operator/(const Grid& a, const Grid& b);
  friend Grid operator/(const Grid& a, double b);
  friend Grid operator/(double a, const Grid& b);
  // *** Printing
  std::string format(char element_delim = ',', char row_delim = '\n', int precision = 4, char end_delim = '\n') const;
  void print() const;
private:
  size_t m_rows;
  size_t m_cols;
  size_t m_size; //m_rows * m_cols
  std::vector<double> m_data;
  Grid ComponentWiseOperation(const Grid& b, const std::function<double(double,double)>& func) const
  {
    assert(m_rows == b.m_rows && m_cols == b.m_cols);
    Grid result(m_rows, m_cols);
    std::vector<double> &result_data = result.m_data;
    const std::vector<double> &b_data = b.m_data;
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0; i<m_size; i++) result_data[i] = func(m_data[i],b_data[i]);
    }
    return result;
  }
  Grid ScalarOperation(double b, const std::function<double(double,double)>& func) const
  {
    Grid result(m_rows, m_cols);
    std::vector<double> &result_data = result.m_data;
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0; i<m_size; i++) result_data[i] = func(m_data[i],b);
    }
    return result;
  }
  Grid UnaryOperation(const std::function<double(double)>& func) const
  {
    Grid result(m_rows, m_cols);
    std::vector<double> &result_data = result.m_data;
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0; i<m_size; i++) result_data[i] = func(m_data[i]);
    }
    return result;
  }
  Grid& InPlaceComponentWiseOperation(const Grid& b, const std::function<double(double&,double)>& func)
  {
    assert(m_rows == b.m_rows && m_cols == b.m_cols);
    const std::vector<double> &b_data = b.m_data;
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0; i<m_size; i++) func(m_data[i],b_data[i]);
    }
    return *this;
  }
  Grid& InPlaceScalarOperation(double b, const std::function<double(double&,double)>& func)
  {
    #pragma omp parallel
    {
      #pragma omp for
      for(int i=0; i<m_size; i++) func(m_data[i],b);
    }
    return *this;
  }
};

#endif