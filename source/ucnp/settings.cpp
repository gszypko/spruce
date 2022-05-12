#include "settings.hpp"

// class constructor
Settings::Settings(fs::path settings_path,std::string unit):
    m_unit_str{unit}   
{
    load_settings(settings_path);
    check_units();
    if (m_runs_found) process_runs();
}

// load .settings file that defines plasma characteristics for MHD simulation
void Settings::load_settings(const fs::path& settings_path)
{
    // mat_in[i] contains [name unit val1 val2 ... valx]
    bool is_settings = settings_path.extension().string()==".settings";
    assert(is_settings && "Path must point to .settings file.");
    std::vector<std::vector<std::string>> mat_in = read_file(settings_path);

    // extract information from input matrix
    m_names.reserve(mat_in.size());  // name for each row
    m_units.reserve(mat_in.size());  // unit for each row
    m_vals.reserve(mat_in.size());   // values associated with each quantity   
    for (auto row : mat_in){
        m_names.push_back(row[0]);
        m_units.push_back(row[1]);
        m_vals.push_back(std::vector(row.begin()+2,row.end()));
    } 

    // determine whether runs are specified
    auto it = std::find(m_names.begin(),m_names.end(),"runs");
    m_runs_found = it != m_names.end();
    if (m_runs_found){
        size_t loc = std::distance(m_names.begin(),it);
        bool single_run_value = m_vals[loc].size() == 1;
        assert(single_run_value && "Only one value can be given for <runs>");
        bool is_positive_int = stoi(m_vals[loc][0]) > 0; 
        assert(is_positive_int && "<runs> must be an integer greater than zero");
        m_runs = stoi(m_vals[loc][0]);
    }
    
    // get unique combination of values for each variable and store in m_unique
    m_unique = unique_comb(m_vals[0],m_vals[1]);
    for (int i = 2; i < m_vals.size(); i++) m_unique = unique_comb(m_unique,m_vals[i]);
}

// verify that all variable units are specified correctly: either <m_unit_str>, <opt>, or another variable name
void Settings::check_units(void) const{
    // possible_units m_unit_str, <opts>, or another variable name
    std::vector<std::string> possible_units(m_names.begin(),m_names.end());
    possible_units.push_back(m_unit_str);
    possible_units.push_back("opt");

    // check that specified units are possible
    for (auto& unit : m_units){
        bool unit_valid = std::find(possible_units.begin(),possible_units.end(),unit) != possible_units.end();
        assert(unit_valid && "Units for each variable must either be <m_unit_str> or another variable name.");
    } 

    // for varialbes with non-m_unit_str units, check that their dependencies are in m_unit_str
    for (auto unit : m_units){
        if (unit!=m_unit_str && unit!="opt"){
            auto it = std::find(m_names.begin(),m_names.end(),unit);
            assert(it!=m_names.end());
            size_t loc = std::distance(m_names.begin(),it);
            bool valid_dependency = m_units[loc]==m_unit_str;
            assert(valid_dependency && "Variables can only be expressed in terms of other variables with m_units_str");
        }
    }    
}

void Settings::process_runs()
{
    // find the variable that corresponds to "runs"
    auto it = std::find(m_names.begin(),m_names.end(),"runs");
    bool runs_found = it!=m_names.end();
    assert(runs_found && "<runs> was not specified within .settings file");
    size_t loc = std::distance(m_names.begin(),it);

    // create unique combinations for each run number
    std::vector<std::vector<std::string>> new_mat;
    new_mat.reserve(m_unique.size()*m_runs);
    for (int i=0; i<m_unique.size(); i++){
        for (int j=0; j<m_runs; j++){
            new_mat.push_back(m_unique[i]);
            new_mat.back()[loc] = num2str(j+1);
        }
    }
    m_unique = new_mat;

}

// choose which set of unique conditions to use
void Settings::choose_array(const int& ind)
{
    assert(ind>=0 && ind<m_unique.size() && "<array> is out of bounds for <m_unique>");
    m_array = m_unique[ind];
    m_array_chosen = true;
}

// return valid task array values
std::string Settings::task_array_range() const
{
    return "Valid task array: [0 "+num2str(array_size()-1)+"]";
}

std::vector<int> Settings::task_array() const
{
    std::vector<int> array(m_unique.size());
    for (int i=0; i<array.size(); i++) array[i] = i;
    return array;
}

// return size of unique array
int Settings::array_size() const
{
    return m_unique.size();
}

// return number of runs
int Settings::runs() const
{
    assert(m_runs_found && "<runs> was not found in the .settings file.");
    return m_runs;
}

// return numeric variable with m_unit_str units
double Settings::getvar(const std::string& name) const
{
    // identify which variable corresponds to <name>
    assert(m_array_chosen && "Set <m_array> using <choose_array> before calling variables with <getvar>.");
    auto it = std::find(m_names.begin(),m_names.end(),name);
    bool found = it != m_names.end();
    assert(found && "Requested name does not correspond to a variable name.");
    size_t loc = std::distance(m_names.begin(),it);
    bool is_opt = m_names[loc] == "opt";
    assert(!is_opt && "Requested variable is of type <opt>, use <getopt> instead.");
    
    // obtain value associated with <name> and convert units if necessary
    double val = stod(m_array[loc]);
    if (m_units[loc]!=m_unit_str){ // if variable units are not "m_unit_str" (i.e., expressed in terms of another variable)
        auto it = std::find(m_names.begin(),m_names.end(),m_units[loc]);
        assert(it!=m_names.end());
        size_t loc2 = std::distance(m_names.begin(),it);
        bool is_numeric = m_units[loc2]==m_unit_str;
        assert(is_numeric && "Variable units can only be expressed in terms of another variable that is expressed in <m_unit_str> units.");
        val *= stod(m_array[loc2]);
    }

    return val;
}

// return option variable as a string
std::string Settings::getopt(const std::string& name) const
{
    assert(m_array_chosen && "Set <m_array> using <choose_array> before calling variables with <getvar>.");
    auto it = std::find(m_names.begin(),m_names.end(),name);
    bool found = it != m_names.end();
    assert(found && "Requested name does not correspond to a variable name");
    size_t loc = std::distance(m_names.begin(),it);
    bool is_opt = m_units[loc] == "opt";
    assert(is_opt && "Requested variable is not of type <opt>");
    return m_array[loc];
}

// write plasma settings for specified m_vals value (i.e., the corresponding row of m_unique)
void Settings::write_array_params(const fs::path& path) const
{
    if (!fs::exists(path)) fs::create_directories(path);
    fs::path file_path{path/"plasma.settings"};
    std::ofstream out_file(file_path);
    for (int i=0; i<m_names.size(); i++){
        out_file << m_names[i] << " = " << m_units[i] << " = " << m_array[i];
        if (i != (m_names.size()-1)) out_file << std::endl;
    }
}

// get unique combinations of each row of mat_in and each element of vec_2
std::vector<std::vector<std::string>> Settings::unique_comb(const std::vector<std::vector<std::string>>& mat_in,const std::vector<std::string>& vec_2) const
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
std::vector<std::vector<std::string>> Settings::unique_comb(const std::vector<std::string>& vec_in,const std::vector<std::string>& vec_2) const
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

// read file with comma delimiter - <#> functions as <,> - all white space ignored - all text in line after <%> is ignored
std::vector<std::vector<std::string>> Settings::read_file(fs::path filePath)
{
    assert(fs::exists(filePath) && !fs::is_empty(filePath) && "File must exist and not be empty.");

    // initialize data with reasonable buffer size
    std::vector<std::vector<std::string>> data;
    data.reserve(100);

    // open input stream to file
    std::ifstream fileStream;       // initialize empty stream object
    while(!fileStream.is_open()){   // while the file is not open...
        fileStream.open(filePath);  // try to open file
    }
    
    // read data from file line-by-line: % comments out line, read after :, spaces are ignored
    while (fileStream.good()){ // while stream is open and has no errors
        std::string currLine;
        std::getline(fileStream,currLine);
        std::size_t pos = currLine.find("%");
        if (pos == 0) currLine.clear();
        else if (pos != std::string::npos) currLine = currLine.substr(0,pos-1);

        // convert '#' to ','
        while ((pos = currLine.find("=")) != std::string::npos) 
            currLine.replace(pos,1,",");

        // remove spaces from string
        auto new_end = std::remove(currLine.begin(),currLine.end(),' ');
        currLine.erase(new_end,currLine.end());
        
        // parse .csv delimited values one at a time
        std::istringstream ss(currLine);
        std::vector<std::string> currVec; // create container for parsing of current line
        currVec.reserve(100); // initialize container size
        while (ss.good()){ // for each delimited value in currLine
            std::string s;
            std::getline(ss,s,','); // place delimited value in 's'
            currVec.push_back(s); // place that element in currVec

            if(currVec.size() == currVec.capacity()) // ensure that reserved size for currVec is large enough
                currVec.reserve(2.*currVec.capacity()); 
        }
        currVec.shrink_to_fit(); // remove excess space

        // store parsed line into data and ensure the reserved storage is large enough
        if (!currLine.empty()){
            data.push_back(currVec); // store in data
            if (data.size() == data.capacity()) 
                data.reserve(2.*data.capacity()); // double buffer size
        }
    }

    // close file stream and shrink data buffer to actual size
    fileStream.close();
    data.shrink_to_fit();

    return data;
}

template <typename T> std::string Settings::num2str(T num,int prec) const
{
    std::stringstream ss;
    ss.precision(prec);
    ss << num;
    return ss.str();
}