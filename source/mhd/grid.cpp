//Class definition for 2D grid of double precision data
//Supports component-wise arithmetic, among other functions
//Built in usage of OpenMP parallelization
#include "grid.hpp"
#define SCALAR_LAMBDA(expr) [](double this_comp, double scalar){return expr;}
#define COMPONENTWISE_LAMBDA(expr) [](double this_comp, double that_comp){return expr;}
#define UNARY_LAMBDA(expr) [=](double this_comp){return expr;}
#define INPLACE_SCALAR_LAMBDA(expr) [](double& this_comp, double scalar){return expr;}
#define INPLACE_COMPONENTWISE_LAMBDA(expr) [](double& this_comp, double that_comp){return expr;}

Grid::Grid(size_t rows, size_t cols): m_rows(rows), m_cols(cols), m_size(rows * cols), m_data(rows * cols) {}
Grid::Grid(size_t rows, size_t cols, double val): m_rows(rows), m_cols(cols), m_size(rows * cols), m_data(rows * cols, val) {}
Grid::Grid(size_t rows, size_t cols, std::vector<double> data): m_rows(rows), m_cols(cols), m_size(rows * cols), m_data(data) {}
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
  assert(curr_min < std::numeric_limits<double>::max() && "Something went wrong in the min calculation");
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
  assert(curr_min < std::numeric_limits<double>::max() && "Something went wrong in the min calculation");
  return curr_min;
}


Grid Grid::min(double b) const { return ScalarOperation( b, SCALAR_LAMBDA(std::min(this_comp,scalar)) ); }

Grid Grid::min(const Grid& b) const { return ComponentWiseOperation( b, COMPONENTWISE_LAMBDA(std::min(this_comp,that_comp)) ); }

double Grid::max() const
{
  double curr_max = -std::numeric_limits<double>::max();
  #pragma omp parallel for reduction(max: curr_max)
  for(int i=0; i<m_size; i++) curr_max = std::max(curr_max,m_data[i]);
  assert(curr_max > -std::numeric_limits<double>::max() && "Something went wrong in the max calculation");
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
  assert(curr_max > -std::numeric_limits<double>::max() && "Something went wrong in the max calculation");
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

Grid Grid::DotProduct2D(const std::vector<Grid>& a, const std::vector<Grid>& b)
{
  assert(a.size() == 2 && b.size() == 2 && "this operator requires 2D vectors");
  return a[0]*b[0] + a[1]*b[1];
}

//Cross product of two vectors in xy-plane (result in z-direction)
Grid Grid::CrossProduct2D(const std::vector<Grid>& a, const std::vector<Grid>& b)
{
  assert(a.size() == 2 && b.size() == 2 && "this operator requires 2D vectors");
  return a[0]*b[1] - a[1]*b[0];
}

//Cross product of vector in z-direction with vector in xy-plane (result in xy-plane)
std::vector<Grid> Grid::CrossProductZ2D(const Grid& a_z, const std::vector<Grid>& b)
{
  assert(b.size() == 2 && "this operator requires 2D vector as second argument");
  Grid b_x = b[0], b_y = b[1];
  return { -a_z*b_y, a_z*b_x };
}

//Cross product of vector in xy-plane with vector in z-direction (result in xy-plane)
std::vector<Grid> Grid::CrossProduct2DZ(const std::vector<Grid>& a, const Grid& b_z)
{
  assert(a.size() == 2 && "this operator requires 2D vector as first argument");
  Grid a_x = a[0], a_y = a[1];
  return Grid::CrossProductZ2D(-1.0*b_z,{a_x,a_y});
}

//General cross product between two vector quantities
std::vector<Grid> Grid::CrossProduct(const std::vector<Grid>& a, const std::vector<Grid>& b)
{
  assert(a.size() == 3 && b.size() == 3 && "this operator requires 3D vectors of Grids");
  Grid a_x = a[0], a_y = a[1], a_z = a[2],
       b_x = b[0], b_y = b[1], b_z = b[2];
  Grid result_x, result_y, result_z;
  result_x = a_y*b_z - a_z*b_y;
  result_y = a_z*b_x - a_x*b_z;
  result_z = a_x*b_y - a_y*b_x;
  return {result_x, result_y, result_z};
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
