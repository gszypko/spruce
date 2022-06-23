#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <iostream>     // std::cerr, std::endl, std::cout
#include <filesystem>   // std::filesystem
#include <algorithm>    // std::remove, std::find, std::distance
#include <vector>       // std::vector
#include <string>       // std::string, std::getline
#include <fstream>      // std::ofstream, std::ifstream
#include <sstream>      // std::istringstream
#include <cassert>      // assert()

namespace fs = std::filesystem;

class Settings
{
public: 
    // *** Alias
    using str_vec = std::vector<std::string>;
    using str_mat = std::vector<std::vector<std::string>>;
    // *** Constructor
    Settings(std::string unit = "cgs"): m_unit_str{unit} {};
    Settings(fs::path settings_path,std::string unit = "cgs");
    static std::unique_ptr<Settings> spawn_settings(std::string name,fs::path settings_path);
    // *** Initialization
    void load_settings(const fs::path& settings_path);
    str_vec get_units_from_vars() const;
    virtual void update_possible_units(){return;}
    void check_units() const;
    void process_runs();
    void choose_array(int array = 0);
    // *** Other Usage
    bool is_unit(const std::string& str) const;
    bool is_name(const std::string& str) const;
    size_t name2ind(const std::string& name) const;
    // *** Getters
    std::vector<int> task_array() const;
    int array_size() const;
    int runs() const;
    str_vec names() const;
    std::string getvar(const std::string& name) const; 
    double getval(const std::string& name) const;
    std::string getopt(const std::string& name) const;
    fs::path set_path(int offset);
    // *** File I/O
    void write_array_params(const fs::path& path,const std::string& name = "file") const;
    std::string task_array_range() const;
    void print_task_array_range() const;
    template <typename T> static std::string num2str(T num,int prec=6)
    {
        std::stringstream ss;
        ss.precision(prec);
        ss << num;
        return ss.str();
    };
protected:
    // *** Members Instantiated at Construction
    const std::string m_unit_str;
    str_vec m_names, m_units;
    str_mat m_vals;
    str_mat m_unique;
    bool m_runs_found; // indicates whether <runs> was given in the .settings 
    int m_runs{-1}; // number of runs for each unique set of conditions
    bool m_array_chosen{false}; // must be set to true by <choose_array> before calling <getvar> or <getopt>
    str_vec m_possible_units;
    // *** Members Handled after Construction
    int m_array{-1};
    str_vec m_array_vals; // m_array_vals[i] is the value correponding to m_names[i] and m_units[i]
    str_mat unique_comb(const str_mat& mat_in,const str_vec& vec_2) const;
    str_mat unique_comb(const str_vec& vec_in,const str_vec& vec_2) const;
    str_mat read_file(fs::path filePath);
};

#endif