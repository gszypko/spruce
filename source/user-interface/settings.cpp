#include "settings.hpp"
#include "ucnp_settings.hpp"

// class constructor
Settings::Settings(fs::path settings_path,std::string unit): m_unit_str{unit}   
{
    load_settings(settings_path);
    if (m_runs_found) process_runs();
    m_possible_units = get_units_from_vars();
    check_units();
    choose_array();
}

// option for constructing a derived class using Settings:: as base
std::unique_ptr<Settings> Settings::spawn_settings(std::string name,fs::path settings_path)
{
    std::unique_ptr<Settings> ptr;
    if (name == "ucnp") 
        ptr = std::unique_ptr<Settings>(new UCNP(settings_path));
    else {
        ptr = std::unique_ptr<Settings>(new Settings(settings_path));
    }
    return ptr;
}

// load .settings file that defines plasma characteristics for MHD simulation
void Settings::load_settings(const fs::path& settings_path)
{
    // mat_in[i] contains [name unit val1 val2 ... valx]
    bool is_settings_file = settings_path.extension().string()==".settings";
    assert(is_settings_file && "Path must point to a file with <.settings> extension.");
    str_mat mat_in = read_file(settings_path);
    // extract information from input matrix
    m_names.reserve(mat_in.size());  // name for each row
    m_units.reserve(mat_in.size());  // unit for each row
    m_vals.reserve(mat_in.size());   // values associated with each quantity   
    for (auto row : mat_in){
        m_names.push_back(row[0]);  // first element is name
        m_units.push_back(row[1]);  // second element is unit
        m_vals.push_back(std::vector(row.begin()+2,row.end())); // all subsequent elements are values associated with (name,unit) pair
    } 
    // determine whether runs is specified and ensure it has a single value that is a positive integer greater than zero
    m_runs_found = is_name("runs");
    if (m_runs_found){
        size_t loc = name2ind("runs");
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

// get units
Settings::str_vec Settings::get_units_from_vars() const
{
    str_vec units = {m_unit_str,"opt"};
    for (int i=0; i<m_names.size(); i++){
        if (m_units[i]==m_unit_str)
            units.push_back(m_names[i]);
    }
    return units;
}

// verify that all variable units are specified correctly: either <m_unit_str>, <opt>, or another variable name
void Settings::check_units(void) const
{
    // for each variable, ensure its unit is one of the possible units
    for (int i=0; i<m_names.size(); i++){
        if (!is_unit(m_units[i])){
            std::cerr << "The units for variable <" + m_names[i] + "> are not valid." << std::endl;
            assert(false);
        }
    }
}

void Settings::process_runs()
{
    assert(is_name("runs"));
    size_t loc = name2ind("runs");
    // create unique combinations for each run number
    str_mat new_mat;
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
void Settings::choose_array(int ind)
{
    assert(ind>=0 && ind<m_unique.size() && "<array> is out of bounds for <m_unique>");
    m_array = ind;
    m_array_vals = m_unique[ind];
    m_array_chosen = true;
}

// check whether <name> is <opt>, m_unit_str, or found within m_names
bool Settings::is_unit(const std::string& str) const
{
    auto it = std::find(m_possible_units.begin(),m_possible_units.end(),str);
    return it != m_possible_units.end();
}

// returns true if <name> is found within <m_names>, otherwise false
bool Settings::is_name(const std::string& str) const
{
    auto it = std::find(m_names.begin(),m_names.end(),str);
    return it != m_names.end();
}

// returns the index of the element of <m_names> that corresponds to <name>
size_t Settings::name2ind(const std::string& name) const
{
    auto it = std::find(m_names.begin(),m_names.end(),name);
    if (it==m_names.end()){
        std::cerr << "Variable name <" << name << "> not found within <m_names>." << std::endl;
        assert(false);
    }
    return std::distance(m_names.begin(),it);
}

void Settings::print_task_array_range() const
{
    std::cout << task_array_range() << std::endl;
}

// return valid task array values
std::string Settings::task_array_range() const
{
    return "Valid task array: [0 "+num2str(array_size()-1)+"]";
}

// returns a vector of valid task arrays
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

Settings::str_vec Settings::names() const {return m_names;}

// return the 
std::string Settings::getvar(const std::string& name) const
{
    assert(m_array_chosen);
    size_t loc,loc2;
    loc = name2ind(name);
    std::string result;
    if (m_units[loc]=="opt") result = getopt(name);
    else {
        result = num2str(getval(name));
    }
    return result;
}

// return option variable as a string
std::string Settings::getopt(const std::string& name) const
{
    assert(m_array_chosen && "Set <m_array_vals> using <choose_array> before calling variables with <getopt>.");
    size_t loc = name2ind(name);
    bool is_opt = m_units[loc] == "opt";
    assert(is_opt && "Requested variable is not of type <opt>");
    return m_array_vals[loc];
}

// return numeric value of the variable associated with <name>
double Settings::getval(const std::string& name) const
{
    // identify which variable corresponds to <name>
    assert(m_array_chosen && "Set <m_array_vals> using <choose_array> before calling variables with <getval>.");
    size_t loc = name2ind(name);
    bool is_opt = m_names[loc] == "opt";
    assert(!is_opt && "Requested variable is of type <opt>, use <getopt> instead.");
    
    // obtain value associated with <name> and convert units if necessary
    double val = stod(m_array_vals[loc]);
    if (m_units[loc]!=m_unit_str){ // if variable units are not "m_unit_str" (i.e., expressed in terms of another variable)
        size_t loc2 = name2ind(m_units[loc]);
        bool is_numeric = m_units[loc2]==m_unit_str;
        assert(is_numeric && "Variable units can only be expressed in terms of another variable that is expressed in <m_unit_str> units.");
        val *= stod(m_array_vals[loc2]);
    }
    return val;
}

fs::path Settings::set_path(int offset)
{
    assert(m_array_chosen);
    int set_num{};
    fs::path result{};
    if (m_runs_found){
        set_num = m_array/m_runs+offset;
        result = "set_" + num2str(set_num);
        result/= "run_" + num2str(getvar("runs"));
    } 
    else{
        set_num = m_array+offset;
        result = "set_" + num2str(set_num);
    } 
    return result;
}

// write plasma settings for specified m_vals value (i.e., the corresponding row of m_unique)
void Settings::write_array_params(const fs::path& path,const std::string& name) const
{
    if (!fs::exists(path)) fs::create_directories(path);
    std::string file_name = name + ".settings";
    fs::path file_path{path/file_name};
    std::ofstream out_file(file_path);
    for (int i=0; i<m_names.size(); i++){
        out_file << m_names[i] << " = " << m_units[i] << " = " << m_array_vals[i];
        if (i != (m_names.size()-1)) out_file << std::endl;
    }
}

// get unique combinations of each row of mat_in and each element of vec_2
Settings::str_mat Settings::unique_comb(const str_mat& mat_in,const str_vec& vec_2) const
{
    str_mat mat_out(mat_in.size()*vec_2.size());
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
Settings::str_mat Settings::unique_comb(const str_vec& vec_in,const str_vec& vec_2) const
{
    str_mat mat_out(vec_in.size()*vec_2.size());
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
Settings::str_mat Settings::read_file(fs::path filePath)
{
    // check that file exists
    assert(fs::exists(filePath) && !fs::is_empty(filePath) && "File must exist and not be empty.");
    // initialize data with reasonable buffer size
    str_mat data;
    data.reserve(100);
    // open input stream to file
    std::ifstream fileStream;       // initialize empty stream object
    while(!fileStream.is_open()){   // while the file is not open...
        fileStream.open(filePath);  // try to open file
    }
    // read data from file line-by-line
    while (fileStream.good()){ // while stream is open and has no errors
        // read in content from current line
        std::string currLine;
        std::getline(fileStream,currLine);
        // remove all content following %
        size_t pos = currLine.find("%");
        if (pos == 0) currLine.clear();
        else if (pos != std::string::npos) currLine = currLine.substr(0,pos-1);
        // convert '#' to ','
        while ((pos = currLine.find("=")) != std::string::npos) 
            currLine.replace(pos,1,",");
        // remove spaces from string
        auto new_end = std::remove(currLine.begin(),currLine.end(),' ');
        currLine.erase(new_end,currLine.end());
        // parse comma-delimited values one at a time
        std::istringstream ss(currLine);
        str_vec currVec; // create container for parsing of current line
        currVec.reserve(100); // initialize container size
        while (ss.good()){
            // get delimited value
            std::string s;
            std::getline(ss,s,',');
            currVec.push_back(s);
            // ensure that currVec has enough capacity to keep pushing back
            if(currVec.size() == currVec.capacity())
                currVec.reserve(2*currVec.capacity()); 
        }
        currVec.shrink_to_fit();
        // store parsed line into data if it is not empty and ensure container size is large enough to keep reading
        if (!currLine.empty()){
            data.push_back(currVec);
            if (data.size() == data.capacity()) 
                data.reserve(2*data.capacity());
        }
    }
    // close file stream and shrink data buffer to actual size
    fileStream.close();
    data.shrink_to_fit();
    return data;
}

