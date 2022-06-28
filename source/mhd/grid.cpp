//Class definition for 2D grid of double precision data
//Supports component-wise arithmetic, among other functions
//Built in usage of OpenMP parallelization
#include "grid.hpp"
#include <limits>
#include <cmath>
#include <algorithm>

#define SCALAR_LAMBDA(expr) [](double this_comp, double scalar){return expr;}
#define COMPONENTWISE_LAMBDA(expr) [](double this_comp, double that_comp){return expr;}
#define UNARY_LAMBDA(expr) [=](double this_comp){return expr;}
#define INPLACE_SCALAR_LAMBDA(expr) [](double& this_comp, double scalar){return expr;}
#define INPLACE_COMPONENTWISE_LAMBDA(expr) [](double& this_comp, double that_comp){return expr;}

Grid::Grid(size_t rows, size_t cols): Grid(rows,cols,0.) {}
Grid::Grid(size_t rows, size_t cols, double val): m_rows(rows), m_cols(cols), m_size(rows * cols), m_data(rows * cols, val) {}
Grid::Grid(size_t rows, size_t cols, std::vector<double> data): m_rows(rows), m_cols(cols), m_size(rows * cols), m_data(data) {}
Grid::Grid(std::vector<double> vec): Grid(1,vec.size(),vec) {}
Grid::Grid(): m_rows(1), m_cols(1), m_size(1), m_data(1) {}
Grid Grid::Zero(size_t rows, size_t cols) { return Grid(rows,cols,0.0); }
Grid Grid::Ones(size_t rows, size_t cols) { return Grid(rows,cols,1.0); }
// generate <num> linearly spaced values between <start> and <end>
Grid Grid::Linspace(double start,double end,int num,int dim)
{
  assert(dim==0 || dim==1);
  assert(end>start);
  assert(num>2);
  Grid result;
  if (dim==0) result = Zero(num,1); // column vector
  if (dim==1) result = Zero(1,num); // row vector
  double spacing = (end-start)/(num-1);
  for (int i=0; i<result.size(); i++) result.m_data[i] = start+i*spacing;
  return result;
}

int Grid::rows() const { return m_rows; }
int Grid::cols() const { return m_cols; }
int Grid::size() const { return m_size; }
std::vector<double> Grid::data() const {return m_data;}
Grid Grid::col(size_t col) const
{
  assert(col>=0 && col<m_cols);
  std::vector<double> vec;
  vec.reserve(m_rows);
  for (int row=0; row<m_rows; row++) vec.push_back((*this)(row,col));
  return Grid(m_rows,1,vec);
  
}

Grid Grid::row(size_t row) const
{
  assert(row>=0 && row<m_rows);
  std::vector<double> vec;
  vec.reserve(m_cols);
  for (int col=0; col<m_cols; col++) vec.push_back((*this)(row,col));
  return Grid(1,m_cols,vec);
}

// find minimum value within the entire grid
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

// returns a grid where all elements smaller than b have been replaced by b
Grid Grid::min(double b) const { return ScalarOperation( b, SCALAR_LAMBDA(std::min(this_comp,scalar)) ); }

// returns a grid where each element is the minimum element
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

Grid Grid::log() const {return UnaryOperation(UNARY_LAMBDA(std::log(this_comp)));}

Grid Grid::sqrt() const { return UnaryOperation( UNARY_LAMBDA(std::sqrt(this_comp)) ); }

Grid Grid::for_each(const Grid& grid,const std::function<double(double,double)>& fun) const {return ComponentWiseOperation(grid,fun);}

Grid Grid::for_each(double a,const std::function<double(double,double)>& fun) const {return ScalarOperation(a,fun);}

// compresses a 2D Grid into a Grid with stacked rows (dim=1) or columns (dim=0)
Grid Grid::vectorise(const Grid& grid,int dim) const
{
  assert(dim==0 || dim==1);
  Grid result;
  std::vector<double> vec;
  vec.reserve(m_size);
  int iter{0};
  if (dim==0){
    for (int i=0; i<m_cols; i++){
      for (int j=0; j<m_rows; j++){
        vec.push_back((*this)(j,i));
      }
    }
    result = Grid(m_size,1,vec);
  } 
  else if (dim==1){
    for (int i=0; i<m_rows; i++){
      for (int j=0; j<m_cols; j++){
        vec.push_back((*this)(i,j));
      }
    }
    result = Grid(1,m_size,vec);
  }
  return result;
}

// 1D integration using the trapezoidal method
Grid Grid::trapz(const Grid& x,const Grid& y)
{
  // determine which dimension is being integrated by inspecting dimensions of <x>
  bool int_row{false},int_col{false};
  if (x.rows()==1) int_row = true;
  else if (x.cols()==1) int_col = true;
  else assert(false);
  // ensure that sizes of <x> and <y> are suitable for integration along that dimension
  int not_integrated,integrated;
  Grid result;
  std::function<double(size_t,size_t)> get_y;
  if (int_row){
    assert(x.cols()==y.cols());
    not_integrated = y.rows();
    integrated = y.cols();
    result = Zero(not_integrated,1);
    get_y = [&](size_t not_integrated,size_t integrated){return y(not_integrated,integrated);};
  }
  if (int_col){
    assert(x.rows()==y.rows());
    not_integrated = y.cols();
    integrated = y.rows();
    result = Zero(1,not_integrated);
    get_y = [&](size_t not_integrated,size_t integrated){return y(integrated,not_integrated);};
  }
  // do integration
  #pragma omp parallel for
  for (int i=0; i<not_integrated; i++){
    for (int j=0; j<integrated-1; j++){
      result(i) += (get_y(i,j)+get_y(i,j+1))*(x(j+1)-x(j))/2.;
    }
  }
  return result;
}

Grid Grid::trapz2D(const Grid& x,const Grid& y,const Grid& z)
{
  assert(x.rows()==z.rows() && x.cols()==1);
  assert(y.cols()==z.cols() && y.rows()==1);
  Grid z_int_x = trapz(x,z);
  Grid z_int_xy = trapz(y,z_int_x);
  return z_int_xy;
}

Grid Grid::trapzcum(const Grid& x, const Grid& y)
{
    // determine which dimension is being integrated by inspecting dimensions of <x>
  bool int_row{false},int_col{false};
  if (x.rows()==1) int_row = true;
  else if (x.cols()==1) int_col = true;
  else assert(false);
  // ensure that sizes of <x> and <y> are suitable for integration along that dimension
  int not_integrated,integrated;
  Grid result = Zero(y.rows(),y.cols());
  std::function<double(size_t,size_t)> get_y;
  if (int_row){
    assert(x.cols()==y.cols());
    not_integrated = y.rows();
    integrated = y.cols();
    get_y = [&](size_t not_integrated,size_t integrated){return y(not_integrated,integrated);};
  }
  if (int_col){
    assert(x.rows()==y.rows());
    not_integrated = y.cols();
    integrated = y.rows();
    get_y = [&](size_t not_integrated,size_t integrated){return y(integrated,not_integrated);};
  }
  // do integration
  #pragma omp parallel for
  for (int i=0; i<not_integrated; i++){
    if (int_row){
      for (int j=0; j<integrated-1; j++){
        result(i,j+1) = result(i,j) + (get_y(i,j)+get_y(i,j+1))*(x(j+1)-x(j))/2.;
      }
    }
    if (int_col){
      for (int j=0; j<integrated-1; j++){
        result(j+1,i) = result(j,i) + (get_y(i,j)+get_y(i,j+1))*(x(j+1)-x(j))/2.;
      }
    }
  }
  return result;
}

Grid Grid::bin_as_list(const Grid& indep,const Grid& dep, const Grid& bin_edges)
{
  assert(bin_edges.size()<dep.size());
  assert(indep.rows()==dep.rows());
  assert(indep.cols()==dep.cols());
  Grid num_in_bin{Zero(1,bin_edges.size()-1)};
  Grid result{Zero(1,bin_edges.size()-1)};
  #pragma omp parallel for
  for (int i=0; i<dep.size(); i++){
    for (int j=0; j<result.size(); j++){
      if ((indep(i)>=bin_edges(j)) && (indep(i)<=bin_edges(j+1))){
        num_in_bin(j)++;
        result(j)+=dep(i);
        break;
      }
    }
  }
  result/=num_in_bin;
  result.m_data.erase(std::remove_if(result.m_data.begin(),result.m_data.end(),[](double x){return !std::isfinite(x);}),result.m_data.end());
  result.m_rows = 1;
  result.m_cols = result.m_data.size();
  result.m_size = result.m_data.size();
  return result;
}

Grid Grid::bin_as_list(const Grid& indep,const Grid& dep, double N, Grid& bins)
{
  bins = Linspace(indep.min(),indep.max(),N);
  double dr = bins(1)-bins(0);
  Grid edges = Linspace(bins.min()-dr/2.,bins.max()+dr/2.,N+1);
  return bin_as_list(indep,dep,edges);
}

Grid Grid::interp_as_list(const Grid& indep,const Grid& dep,const Grid& query)
{
  assert(indep.size()==dep.size());
  for (auto& val : query.m_data) assert(val>=indep.min() && val<=indep.max());
  Grid result = Zero(query.rows(),query.cols());
  #pragma omp parallel for
  for (int i=0; i<query.size(); i++){
    double curr_min = -std::numeric_limits<double>::max();
    double curr_max = std::numeric_limits<double>::max();
    std::vector<double> min_vals{}, max_vals{}, equal_vals{};
    equal_vals.reserve(indep.size());
    min_vals.reserve(indep.size());
    max_vals.reserve(indep.size());
    if (indep(0)==query(i)) equal_vals.push_back(dep(0));
    for (int j=0; j<indep.size(); j++){
      if (indep(j)==query(i)) equal_vals.push_back(dep(j));
      else if (indep(j)==curr_min) min_vals.push_back(dep(j));
      else if (indep(j)==curr_max) max_vals.push_back(dep(j));
      else if (indep(j)<query(i) && indep(j)>curr_min){
        curr_min = indep(j);
        min_vals.clear();
        min_vals.push_back(dep(j));
      }
      else if (indep(j)>query(i) && indep(j)<curr_max){
        curr_max = indep(j);
        max_vals.clear();
        max_vals.push_back(dep(j));
      }
    }
    if (equal_vals.size()>0){
      for (auto& val : equal_vals) result(i)+=val/equal_vals.size();
    }
    else{
      double min{}, max{};
      for (auto& val : min_vals) min+=val/min_vals.size();
      for (auto& val : max_vals) max+=val/max_vals.size();
      result(i) = min + (query(i)-curr_min)*(max-min)/(curr_max-curr_min);
    }
  }
  return result;
}

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

void Grid::print() const
{
  for (int i=0; i<m_rows; i++){
    std::cout << "  ";
    for (int j=0; j<m_cols; j++){
      std::cout << (*this)(i,j) << " ";
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
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

double& Grid::operator()(size_t el)
{
  assert(el>=0 && el<m_size);
  return m_data[el];
}

double Grid::operator()(size_t el) const
{
  assert(el>=0 && el<m_size);
  return m_data[el];
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
