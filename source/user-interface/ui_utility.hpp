#ifndef UI_UTILITY_HPP
#define UI_UTILITY_HPP

#include "constants.hpp" // global macro definitions of constants in cgs units
#include "grid.hpp" // Grid:: parallelized matrix class used in MHD simulation
#include "PlasmaSettings.hpp" // class definition for plasma definition
#include "MhdInp.hpp" // input parameters for MHD solver

#include <iostream> // std::cout, std::endl
#include <filesystem> // std::filesystem - for file directory interactivity
#include <fstream> // std::ofstream, std::ifstream - file i/o
#include <sstream> // std::istringstream
#include <vector> // std::vector
#include <string> // std::string, std::getline
#include <math.h> // sqrt, pow, etc - basic math functions
#include <algorithm> // std::remove

//*** GENERAL FUNCTIONS ***//

int round2int(const double& a);
bool strcmp(std::string str1, std::string str2);
std::vector<double> bin_vector(const std::vector<double>& vec_in,const std::vector<double>& bins);
double euclidean_norm(const std::vector<double>& vec_in);
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::vector<std::string>>& mat_in,const std::vector<std::string>& vec_2);
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::string>& vec_in,const std::vector<std::string>& vec_2);
std::vector<std::vector<std::string>> readCSV(std::filesystem::path filePath);

//*** FUNCTIONS FOR USE WITH GRID.HPP ***//

void meshgrid(const std::vector<double>& v_x,const std::vector<double>& v_y,Grid& grid_x,Grid& grid_y);
Grid gaussian2D(Grid x,Grid y,double amp,double sig_x,double sig_y,double x_cen,double y_cen);
Grid exponential2D(Grid x,Grid y,double amp,double sig_x,double sig_y,double x_cen,double y_cen);
std::vector<double> quadrupole_field(double x,double y,double z,double Bx);

//*** FUNCTION TEMPLATES ***//

// generate linearly spaced values between min and max
template <typename T> std::vector<T> linspace(double min, double max, int num_pts)
{
    std::vector<T> vec_out(num_pts);
    if (num_pts < 1) std::cerr << "Error: must request more than one bin with linspace";
    double spacing = (max - min)/(num_pts-1);
    for (int i = 0; i < num_pts; i++) vec_out[i] = min + spacing*i;
    return vec_out;
}

// convert a double or integer to string with specified precision
template <typename T> std::string num2str(T num,int prec = 3)
{
    std::stringstream ss;
    ss.precision(prec);
    ss << num;

    return ss.str();
}

// append row to csv
template <typename T> void append_row_to_csv(const std::filesystem::path file_path, const std::vector<T> data)
// filePath: full path to .csv file specified as type std::filesystem::path
// data: a std::vector of fundamental data type (i.e., string, double, int, etc)
{
    // open output stream to file depending on whether file exists or not
    std::ofstream out_file;
    if (!std::filesystem::exists(file_path)){ // if file does not exist
        out_file.open(file_path); // open stream to file
    }
    else if (std::filesystem::exists(file_path)){ // if file does exist
        out_file.open(file_path,std::ofstream::app);
    }

    // append 'data' to file
    std::string delim = ",";
    for (int i = 0; i < data.size() - 1; i++){ // for each element, except for the last one
        out_file << data[i] << delim;
    }
    out_file << data.back() << std::endl; // use end line for the last element so that future writing automatically goes to new line
}

//*** PHYSICS FORMULAS - all use cgs units ***//

double coulomb_coupling(const double& n,const double& T);
double wigner_seitz_radius(const double& n);
double plasma_freq(const double& n,const double& m);
double einstein_freq(const double& n,const double& m);
double nearest_coulomb_pot(const double& n);
double debye_length(const double& n,const double& T);
double screening_parameter(const double& n, const double& Te);
double dih_temp(const double& n,const double& Te);
double get_tau_exp(const double& sig,const double& m, const double& Te);

#endif