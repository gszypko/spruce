#include "ui_utility.hpp"

// round double to int
int round2int(const double& a) { return static_cast<int>(a < 0 ? a - 0.5 : a + 0.5); }

// bin <vec_in> into <num_bins> linearly spaced bins between min and max
std::vector<double> bin_vector(const std::vector<double>& vec_in,const std::vector<double>& bins)
{
    std::vector<double> num_in_bin(bins.size(),0.);
    double bin_spacing = bins[1] - bins[0];
    
    #pragma omp parallel for
    for (int i = 0; i < bins.size(); i++){ 
        for (int j = 0; j < vec_in.size(); j++){
            bool cond1, cond2, in_bin;
            cond1 = vec_in[j] < (bins[i] + bin_spacing/2.); 
            cond2 = vec_in[j] > (bins[i] - bin_spacing/2.); 
            in_bin = cond1 && cond2;
            if (in_bin) num_in_bin[i] += 1./vec_in.size();
        }
    }
    return num_in_bin;
}

// get unique combinations of each row of mat_in and each element of vec_2
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::vector<std::string>>& mat_in,const std::vector<std::string>& vec_2)
{
    std::vector<std::vector<std::string>> mat_out(mat_in.size()*vec_2.size());
    int row = 0;
    for (int i = 0; i < mat_in.size(); i++){
        for (int j = 0; j < vec_2.size(); j++){
            mat_out[row].reserve(mat_in[i].size()+1);
            for (int k = 0; k < mat_in[i].size(); k++) mat_out[row].push_back(mat_in[i][k]);
            mat_out[row].push_back(vec_2[j]);
            row++;
        }
    }

    return mat_out;
}

// get unique combinations of the elements in the two input vectors
std::vector<std::vector<std::string>> unique_comb(const std::vector<std::string>& vec_in,const std::vector<std::string>& vec_2)
{
    std::vector<std::vector<std::string>> mat_out(vec_in.size()*vec_2.size());
    int row = 0;
    for (int i = 0; i < vec_in.size(); i++){
        for (int j = 0; j < vec_2.size(); j++){
            mat_out[row].push_back(vec_in[i]);
            mat_out[row].push_back(vec_2[j]);
            row++;
        }
    }

    return mat_out;
}

// compute the norm (i.e., square root of a vectors inner product with itself)
double euclidean_norm(const std::vector<double>& vec_in)
{
    double norm;
    for (auto val : vec_in) norm += pow(val,2.);
    norm = sqrt(norm);
    return norm;
}

// read file with comma delimiter, "=" functions as ",", all text after "%" is ignored
std::vector<std::vector<std::string>> readCSV(std::filesystem::path filePath)
// filePath: full path to .csv file to read from, formatted as std::filesystem::path
{
    // check that file exists and is not empty
    if (!std::filesystem::exists(filePath) || std::filesystem::is_empty(filePath)){
        std::cerr << "The specified .csv file is either empty or does not exist." << std::endl;
    }

    // initialize data with reasonable buffer size
    std::vector<std::vector<std::string>> data;
    data.reserve(100);

    // open input stream to file
    std::ifstream fileStream;            // initialize empty stream object
    while(!fileStream.is_open()){             // while the file is not open...
        fileStream.open(filePath);  // try to open file
    }
    
    // read data from file line-by-line: % comments out line, read after :, spaces are ignored
    while (fileStream.good()){ // while stream is open and has no errors
        // read in current line from .csv file
        std::string currLine; // initialize container for current line of .csv file
        std::getline(fileStream,currLine); // read in current line

        // if first character is "%", line is ignored
        std::size_t pos = currLine.find("%");
        if (pos == 0) currLine.clear();
        else if (pos != std::string::npos) currLine = currLine.substr(0,pos-1);

        // convert '#' to ','
        while ((pos = currLine.find("=")) != std::string::npos) currLine.replace(pos,1,",");

        // remove spaces from string
        currLine.erase(std::remove(currLine.begin(),currLine.end(),' '),currLine.end());
        
        std::istringstream ss(currLine); // create stream to current line

        // parse .csv delimited values one at a time
        std::vector<std::string> currVec; // create container for parsing of current line
        currVec.reserve(100); // initialize container size
        while (ss.good()){ // for each delimited value in currLine
            std::string s;
            std::getline(ss,s,','); // place delimited value in 's'
            currVec.push_back(s); // place that element in currVec

            if(currVec.size() == currVec.capacity()){ // ensure that reserved size for currVec is large enough
                currVec.reserve(2.*currVec.capacity()); 
            }
        }
        currVec.shrink_to_fit(); // remove excess space

        // store parsed line into data and ensure the reserved storage is large enough
        
        if (!currLine.empty()){
            data.push_back(currVec); // store in data
            if(data.size() == data.capacity()){ // if actual size has reached buffer size
            data.reserve(2.*data.capacity()); // double buffer size
        }
        
    }
    }

    // close file stream and shrink data buffer to actual size
    fileStream.close();
    data.shrink_to_fit();

    return data;
}

// compare <str1> and <str2>
bool strcmp(std::string str1, std::string str2)
{
    bool is_equal;
    if (str1.compare(str2) == 0) is_equal = true;
    else is_equal = false;
    return is_equal;
}

// return <x> and <y> mesh grids 
void meshgrid(const std::vector<double>& v_x,const std::vector<double>& v_y,Grid& grid_x,Grid& grid_y)
{
    grid_x = Grid(v_x.size(),v_y.size(),0.);
    grid_y = Grid(v_x.size(),v_y.size(),0.);

    for (int i = 0; i < v_x.size(); i++){
        for (int j = 0; j < v_y.size(); j++){
            grid_x(i,j) = v_x[i];
            grid_y(i,j) = v_y[j];
        }
    }
}

// define 2D gaussian distribution - spatial units must match input grid, output grid units match <amp>
Grid gaussian2D(Grid x,Grid y,double amp,double min,double sig_x,double sig_y,double x_cen,double y_cen)
// x: input meshgrid of x-positions
// y: input meshgrid of y-positions
// amp: amplitude of 2D Gaussian, can be +/-
// sig_x: RMS width of Gaussian distribution on x-axis
// sig_y: RMS width of Gaussian distribution on y-axis
// x_cen: center of Gaussian distribution on x-axis
// y_cen: center of Gaussian distribution on y-axis
// All position units must be the same, otherwise they don't matter. Output grid has same units as <amp>.
{
    bool size_check { true };
    if (x.rows() != y.rows()) size_check = false;
    if (x.cols() != y.cols()) size_check = false;
    if (x.size() != y.size()) size_check = false;
    if (!size_check) std::cerr << "Input grids <x> and <y> must have the same size." << std::endl;
    Grid g2D { Grid(x.rows(),x.cols()) };
    for (int i = 0; i < g2D.rows(); i++){
        for (int j = 0; j < g2D.cols(); j++){
            g2D(i,j) = min+amp*exp(-0.5*pow((x(i,j)-x_cen)/abs(sig_x),2))*exp(-0.5*pow((y(i,j)-y_cen)/abs(sig_y),2));
        }
    }
    return g2D;
}

// define 2D exponential distribution - spatial units must match input grid, output grid units match <amp>
Grid exponential2D(Grid x,Grid y,double amp,double min,double sig_x,double sig_y,double x_cen,double y_cen)
// x: input meshgrid of x-positions
// y: input meshgrid of y-positions
// amp: amplitude of 2D Gaussian, can be +/-
// sig_x: RMS width of Gaussian distribution on x-axis
// sig_y: RMS width of Gaussian distribution on y-axis
// x_cen: center of Gaussian distribution on x-axis
// y_cen: center of Gaussian distribution on y-axis
// All position units must be the same, otherwise they don't matter. Output grid has same units as <amp>.
{
    bool size_check { true };
    if (x.rows() != y.rows()) size_check = false;
    if (x.cols() != y.cols()) size_check = false;
    if (x.size() != y.size()) size_check = false;
    if (!size_check) std::cerr << "Input grids <x> and <y> must have the same size." << std::endl;
    Grid exp2D { Grid(x.rows(),x.cols()) };
    for (int i = 0; i < exp2D.rows(); i++){
        for (int j = 0; j < exp2D.cols(); j++){
            exp2D(i,j) = min+amp*exp(-abs(x(i,j)-x_cen)/abs(sig_x))*exp(-abs(y(i,j)-y_cen)/abs(sig_y));
        }
    }
    return exp2D;
}

// get quadrupole magnetic field, cgs units, symmetry axis = x
std::vector<double> quadrupole_field(double x,double y,double z,double dBdx)
{
    enum Ax { ax_x, ax_y, ax_z , ax_num};
    std::vector<double> B(ax_num);
    B[ax_x] = -dBdx*x;
    B[ax_y] = dBdx*y/2.;
    B[ax_z] = dBdx*z/2.;
    return B;
}

// ratio between average Coulomb interaction energy and thermal energy (inputs in cgs)
double coulomb_coupling(const double& n,const double& T)
// n: plasma density (cm^-3)
// T: temperature (K) of a particular species (i.e., electron or ion)
{ return nearest_coulomb_pot(n)/(K_B*T); }

// average interparticle spacing of gas with density n (SI units)
double wigner_seitz_radius(const double& n) 
// n: plasma density (cm^-3)
{ return pow(3./(4.*PI*n),1./3.); } 

// plasma oscillation frequency (rad/s, inputs cgs)
double plasma_freq(const double& n,const double& m) 
// n: plasma density (cm^-3)
// m: particle mass (g)
{ return sqrt(4*PI*n*pow(E,2.)/m); }

// einsten frequency (rad/s, inputs cgs)
double einstein_freq(const double& n,const double& m)
// n: plasma density (cm^-3)
// m: particle mass (g)
{ return plasma_freq(n,m)/sqrt(3.); }

// average Coulomb potential between nearest neighbors, cgs units
double nearest_coulomb_pot(const double& n) 
// n: plasma density (cm^-3)
{ return E*E/wigner_seitz_radius(n); } 

// debye length of a particular plasma species, cgs units
double debye_length(const double& n,const double& T) 
// n: plasma density (cm^-3)
// T: temperature (K) of a particular plasma species
{ return sqrt(K_B*T/(4*PI*n*E*E));}

// plasma screening parameter for electrons, dimensionless, inputs cgs
double screening_parameter(const double& n, const double& Te)
// n: plasma density (cm^-3)
// Te: electron temperature (K)
{ return wigner_seitz_radius(n)/debye_length(n,Te); } 

// equilibrium plasma temperature due to disorder-induced heating, inputs cgs
double dih_temp(const double& n,const double& Te) 
// n: plasma density (cm^-3)
// Te: electron temperature (K)
{ return 2.*nearest_coulomb_pot(n)/3./K_B*(1.+screening_parameter(n,Te)/2.); }

// timescale for a Gaussian UCNP expansion into vacuum, cgs units
double tau_exp(const double& sig,const double& m, const double& Te)
// n: plasma density (cm^-3)
// m: ion mass (g)
// Te: electron temperature (K)
{
    return sqrt(m*pow(sig,2.)/(K_B*Te));
}