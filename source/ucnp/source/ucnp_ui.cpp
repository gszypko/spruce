#include "ucnp_ui.hpp"

// round double to int
int round2int(const double a) { return static_cast<int>(a < 0 ? a - 0.5 : a + 0.5); }

// bin vec_in into num_bins between min and max
std::vector<double> bin_vector(std::vector<double> vec_in, std::vector<double> bins)
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
std::vector<std::vector<double>> unique_comb(const std::vector<std::vector<double>>& mat_in,const std::vector<double>& vec_2)
{
    std::vector<std::vector<double>> mat_out(mat_in.size()*vec_2.size());
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
std::vector<std::vector<double>> unique_comb(const std::vector<double>& vec_in,const std::vector<double>& vec_2)
{
    std::vector<std::vector<double>> mat_out(vec_in.size()*vec_2.size());
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
double euclidean_norm(std::vector<double> vec_in)
{
    double norm;
    for (auto val : vec_in) norm += pow(val,2.);
    norm = sqrt(norm);
    return norm;
}

namespace fs = std::filesystem;
using namespace std;

// read file with comma delimiter, % used as comment entire line
std::vector<std::vector<std::string>> readCSV(std::filesystem::path filePath)
// filePath: full path to .csv file to read from, formatted as std::filesystem::path
{
    // check that file exists and is not empty
    if (!std::filesystem::exists(filePath) || std::filesystem::is_empty(filePath)){
        cerr << "The specified .csv file is either empty or does not exist." << endl;
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
        else if (pos != string::npos) currLine = currLine.substr(0,pos-1);

        // remove spaces from string
        // currLine.erase(std::remove(currLine.begin(),currLine.end(),' '),currLine.end());

        // all text before ":" on any line is ignored
        pos = currLine.find("=");
        if (pos != string::npos) currLine = currLine.substr(pos+1,string::npos);
        
        std::istringstream ss(currLine); // create stream to current line

        // parse .csv delimited values one at a time
        std::vector<std::string> currVec; // create container for parsing of current line
        currVec.reserve(100); // initialize container size
        while (ss.good()){ // for each delimited value in currLine
            string s;
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

// ratio between average Coulomb interaction energy and thermal energy (dimensionless)
double coulomb_coupling(const double n,const double T)
// n: plasma density (m^-3)
// T: temperature (K) of a particular species (i.e., electron or ion)
{ return cts::e*cts::e/(4.*cts::pi*cts::eps0*wigner_seitz_radius(n)*cts::kB*T); }

// average interparticle spacing of gas with density n (SI units)
double wigner_seitz_radius(const double n) 
// n: plasma density (m^-3)
{ return pow(3./(4.*cts::pi*n),1./3.); } 

// plasma oscillation frequency (rad/s)
double plasma_freq(const double n,const double m) 
// n: plasma density (m^-3)
// m: particle mass (kg)
{ return sqrt(n*pow(cts::e,2.)/(m*cts::eps0)); }

// einsten frequency (rad/s)
double einstein_freq(const double n,const double m)
// n: plasma density (m^-3)
// m: particle mass (kg)
{ return plasma_freq(n,m)/sqrt(3.); }

// average Coulomb potential between nearest neighbors (j)
double nearest_coulomb_pot(const double n) 
// n: plasma density (m^-3)
{ return cts::e*cts::e/(4.*cts::pi*cts::eps0*wigner_seitz_radius(n)); } 

// debye length of a particular plasma species
double debye_length(const double n,const double T) 
// n: plasma density (m^-3)
// T: temperature (K) of a particular plasma species
{ return sqrt(cts::eps0*cts::kB*T/(n*cts::e*cts::e));}

// plasma screening parameter for electrons
double screening_parameter(const double n, const double Te)
// n: plasma density (m^-3)
// Te: electron temperature (K)
{ return wigner_seitz_radius(n)/debye_length(n,Te); } 

// equilibrium plasma temperature due to disorder-induced heating
double dih_temp(const double n,const double Te) 
// n: plasma density (m^-3)
// Te: electron temperature (K)
{ return 2.*nearest_coulomb_pot(n)/3./cts::kB*(1.+screening_parameter(n,Te)/2.); }