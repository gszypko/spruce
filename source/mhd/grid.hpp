//Class definition for 2D grid of double precision data
//Supports component-wise arithmetic, among other functions
//Built in usage of OpenMP parallelization
#ifndef GRID_HPP
#define GRID_HPP

#include "instrumentor.hpp"
#include <vector>
#include <string>
#include <functional>
#include <cassert>
#include <omp.h>
#include <sstream>
#include <limits>
#include <cmath>

class Grid
{
public:
  Grid(size_t rows, size_t cols);
  Grid(size_t rows, size_t cols, double val);
  Grid();
  static Grid Zero(size_t rows, size_t cols);
  static Grid Ones(size_t rows, size_t cols);
  double rows() const;
  double cols() const;
  double size() const;
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
  std::string format(char element_delim = ',', char row_delim = '\n', int precision = 4, char end_delim = '\n') const;
  static std::vector<Grid> CrossProduct(const Grid& a_x, const Grid& a_y, const Grid& a_z, const Grid& b_x, const Grid& b_y, const Grid& b_z);
  // OPERATOR OVERRIDES
  double& operator()(size_t i, size_t j);
  double operator()(size_t i, size_t j) const;
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
      #if BENCHMARKING_ON
      InstrumentationTimer timer("componentwise op loop thread");
      #endif
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
      #if BENCHMARKING_ON
      InstrumentationTimer timer("scalar op loop thread");
      #endif
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
      #if BENCHMARKING_ON
      InstrumentationTimer timer("unary op loop thread");
      #endif
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
      #if BENCHMARKING_ON
      InstrumentationTimer timer("in place componentwise op loop thread");
      #endif
      #pragma omp for
      for(int i=0; i<m_size; i++) func(m_data[i],b_data[i]);
    }
    return *this;
  }
  Grid& InPlaceScalarOperation(double b, const std::function<double(double&,double)>& func)
  {
    #pragma omp parallel
    {
      #if BENCHMARKING_ON
      InstrumentationTimer timer("in place scalar op loop thread");
      #endif
      #pragma omp for
      for(int i=0; i<m_size; i++) func(m_data[i],b);
    }
    return *this;
  }
};

#endif