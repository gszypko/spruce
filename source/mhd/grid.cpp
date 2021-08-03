//Class definition for 2D grid of double precision data
//Supports component-wise arithmetic, among other functions
//Built in usage of OpenMP parallelization
#include <vector>
#include <assert.h>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include "grid.hpp"
#define SCALAR_LAMBDA(expr) [](double this_comp, double scalar){return expr;}
#define COMPONENTWISE_LAMBDA(expr) [](double this_comp, double that_comp){return expr;}
#define UNARY_LAMBDA(expr) [=](double this_comp){return expr;}
#define INPLACE_SCALAR_LAMBDA(expr) [](double& this_comp, double scalar){return expr;}
#define INPLACE_COMPONENTWISE_LAMBDA(expr) [](double& this_comp, double that_comp){return expr;}

Grid::Grid(size_t rows, size_t cols): m_rows(rows), m_cols(cols), m_size(rows * cols), m_data(rows * cols) {}
Grid::Grid(size_t rows, size_t cols, double val): m_rows(rows), m_cols(cols), m_size(rows * cols), m_data(rows * cols, val) {}
Grid::Grid(): m_rows(1), m_cols(1), m_size(1), m_data(1) {}
Grid Grid::Zero(size_t rows, size_t cols) { return Grid(rows,cols,0.0); }
Grid Grid::Ones(size_t rows, size_t cols) { return Grid(rows,cols,1.0); }

double Grid::rows() const { return m_rows; }
double Grid::cols() const { return m_cols; }
double Grid::size() const { return m_size; }

double Grid::min() const
{
  double curr_min = std::numeric_limits<double>::max();
  #pragma omp parallel for reduction(min: curr_min)
  for(int i=0; i<m_size; i++) curr_min = std::min(curr_min,m_data[i]);
  return curr_min;
}

//Finds the min within the rectangular range defined by the bounds (il,jl) and (iu,ju)
double Grid::min(int il, int jl, int iu, int ju) const
{
  double curr_min = std::numeric_limits<double>::max();
  #pragma omp parallel for collapse(2) reduction(min: curr_min)
  for(int i=il; i<=iu; i++){
    for(int j=jl; j<=ju; j++){
      curr_min = std::min(curr_min,(*this)(i,j));
    }
  }
  return curr_min;
}


Grid Grid::min(double b) const { return ScalarOperation( b, SCALAR_LAMBDA(std::min(this_comp,scalar)) ); }

Grid Grid::min(const Grid& b) const { return ComponentWiseOperation( b, COMPONENTWISE_LAMBDA(std::min(this_comp,that_comp)) ); }

double Grid::max() const
{
  double curr_max = -std::numeric_limits<double>::max();
  #pragma omp parallel for reduction(max: curr_max)
  for(int i=0; i<m_size; i++) curr_max = std::max(curr_max,m_data[i]);
  return curr_max;
}

//Finds the min within the rectangular range defined by the bounds (il,jl) and (iu,ju)
double Grid::max(int il, int jl, int iu, int ju) const
{
  double curr_max = -std::numeric_limits<double>::max();
  #pragma omp parallel for collapse(2) reduction(max: curr_max)
  for(int i=il; i<=iu; i++){
    for(int j=jl; j<=ju; j++){
      curr_max = std::max(curr_max,(*this)(i,j));
    }
  }
  return curr_max;
}


Grid Grid::max(double b) const { return ScalarOperation( b, SCALAR_LAMBDA(std::max(this_comp,scalar)) ); }

Grid Grid::max(const Grid& b) const { return ComponentWiseOperation( b, COMPONENTWISE_LAMBDA(std::max(this_comp,that_comp)) ); }

Grid Grid::square() const { return (*this)*(*this); }

Grid Grid::abs() const { return UnaryOperation( UNARY_LAMBDA(std::abs(this_comp)) ); }

Grid Grid::pow(double power) const { return UnaryOperation( UNARY_LAMBDA(std::pow(this_comp, power)) ); }

Grid Grid::sqrt() const { return UnaryOperation( UNARY_LAMBDA(std::sqrt(this_comp)) ); }

//Formats the grid as a character-delimited string, with num of sig figs given by precision
//Default precision is 4 sig figs; setting precision to -1 produces maximum (stable) precision
std::string Grid::format(char element_delim, char row_delim, int precision, char end_delim) const
{
  std::ostringstream out;
  if(precision == -1) out.precision(std::numeric_limits<double>::digits10 + 1);
  else out.precision(precision);
  for(int i=0; i<m_rows; i++){
    for(int j=0; j<m_cols; j++){
      out << (*this)(i,j) << element_delim;
    }
    out.seekp(-1, std::ios_base::end);
    out << row_delim;
  }
  out.seekp(-1, std::ios_base::end);
  out << end_delim;
  return out.str();
}

double& Grid::operator()(size_t i, size_t j)
{
  assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
  return m_data[i*m_cols + j];
}

double Grid::operator()(size_t i, size_t j) const
{
  assert(i >= 0 && i < m_rows && j >= 0 && j < m_cols);
  return m_data[i*m_cols + j];
}

Grid& Grid::operator+=(const Grid& b) { return (*this).InPlaceComponentWiseOperation( b, INPLACE_COMPONENTWISE_LAMBDA(this_comp += that_comp) ); }

Grid& Grid::operator+=(double b) { return (*this).InPlaceScalarOperation( b, INPLACE_SCALAR_LAMBDA(this_comp += scalar) ); }

Grid& Grid::operator-=(const Grid& b) { return (*this).InPlaceComponentWiseOperation( b, INPLACE_COMPONENTWISE_LAMBDA(this_comp -= that_comp) ); }

Grid& Grid::operator-=(double b) { return (*this).InPlaceScalarOperation( b, INPLACE_SCALAR_LAMBDA(this_comp -= scalar) ); }

Grid& Grid::operator*=(const Grid& b) { return (*this).InPlaceComponentWiseOperation( b, INPLACE_COMPONENTWISE_LAMBDA(this_comp *= that_comp) ); }

Grid& Grid::operator*=(double b) { return (*this).InPlaceScalarOperation( b, INPLACE_SCALAR_LAMBDA(this_comp *= scalar) ); }

Grid& Grid::operator/=(const Grid& b) { return (*this).InPlaceComponentWiseOperation( b, INPLACE_COMPONENTWISE_LAMBDA(this_comp /= that_comp) ); }

Grid& Grid::operator/=(double b) { return (*this).InPlaceScalarOperation( b, INPLACE_SCALAR_LAMBDA(this_comp /= scalar) ); }

std::ostream& operator<<(std::ostream& os, const Grid& a) { os << a.format(); return os; }

Grid operator+(const Grid& a, const Grid& b) { return a.ComponentWiseOperation( b, COMPONENTWISE_LAMBDA(this_comp + that_comp) ); }

Grid operator+(const Grid& a, double b) { return a.ScalarOperation( b, SCALAR_LAMBDA(this_comp + scalar) ); }

Grid operator+(double a, const Grid& b) { return b+a; }

Grid operator-(const Grid& a, const Grid& b) { return a.ComponentWiseOperation( b, COMPONENTWISE_LAMBDA(this_comp - that_comp) ); }

Grid operator-(const Grid& a, double b) { return a.ScalarOperation( b, SCALAR_LAMBDA(this_comp - scalar) ); }

Grid operator-(double a, const Grid& b) { return b.ScalarOperation( a, SCALAR_LAMBDA(scalar - this_comp) ); }

Grid operator-(const Grid& a) { return -1.0*a; }

Grid operator*(const Grid& a, const Grid& b) { return a.ComponentWiseOperation( b, COMPONENTWISE_LAMBDA(this_comp * that_comp) ); }

Grid operator*(const Grid& a, double b) { return a.ScalarOperation( b, SCALAR_LAMBDA(this_comp * scalar) ); }

Grid operator*(double a, const Grid& b) { return b*a; }

Grid operator/(const Grid& a, const Grid& b) { return a.ComponentWiseOperation( b, COMPONENTWISE_LAMBDA(this_comp / that_comp) ); }

Grid operator/(const Grid& a, double b) { return a.ScalarOperation( b, SCALAR_LAMBDA(this_comp / scalar) ); }

Grid operator/(double a, const Grid& b) { return b.ScalarOperation( a, SCALAR_LAMBDA(scalar / this_comp) ); }
