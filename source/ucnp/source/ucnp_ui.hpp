#ifndef UCNPUI_HPP
#define UCNPUI_HPP

#include "grid.hpp"
#include "plasmadomain.hpp"

#include <iostream>
#include <vector>
#include <string>                   // std::string
#include <algorithm>                // std::remove
#include <math.h>
#include <filesystem>
#include <sstream>
#include <constants_ucnp.hpp>

//*** General Functions ***//

int round2int(const double a);
std::vector<double> bin_vector(std::vector<double> vec_in,std::vector<double> bins);
double euclidean_norm(std::vector<double> vec_in);
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::vector<std::string>>& mat_in,const std::vector<std::string>& vec_2);
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::string>& vec_in,const std::vector<std::string>& vec_2);
std::vector<std::vector<std::string>> readCSV(std::filesystem::path filePath);
//*** Physics Formulas ***//

double coulomb_coupling(const double n,const double T);
double wigner_seitz_radius(const double n);
double plasma_freq(const double n,const double m);
double einstein_freq(const double n,const double m);
double nearest_coulomb_pot(const double n);
double debye_length(const double n,const double T);
double screening_parameter(const double n, const double Te);
double dih_temp(const double n,const double Te);

//*** Templates ***//

// generate linearly spaced values between min and max
template <typename T> 
std::vector<T> linspace(double min, double max, int num_pts)
{
    std::vector<T> vec_out(num_pts);
    if (num_pts < 1) std::cerr << "Error: must request more than one bin with linspace";
    double spacing = (max - min)/(num_pts-1);
    for (int i = 0; i < num_pts; i++) vec_out[i] = min + spacing*i;
    return vec_out;
}

// convert a double or integer to string with specified precision
template <typename T>
std::string num2str(T num,int prec = 3)
{
    std::stringstream ss;
    ss.precision(prec);
    ss << num;

    return ss.str();
}

// append row to csv
template <typename T>
void append_row_to_csv(const std::filesystem::path file_path, const std::vector<T> data)
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


#endif