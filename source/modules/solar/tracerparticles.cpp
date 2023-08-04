#include "tracerparticles.hpp"
#include "plasmadomain.hpp"
#include "solarutils.hpp"
#include "idealmhd.hpp"
#include "utils.hpp"
#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>

TracerParticles::TracerParticles(PlasmaDomain &pd) : Module(pd) {
  assert(dynamic_cast<IdealMHD*>(m_pd.m_eqs.get()) && \
    "Module designed for IdealMHD EquationSet (ensure that equation_set is set before modules in the config)");
}

void TracerParticles::setupModule(){
  assert(m_particles.size() <= m_pd.m_xdim && "You cannot have more than xdim particles in a single TracerParticles module");
  for(int i=0; i<m_pd.m_xdim; i++) x_vec.push_back(m_pd.grid(PlasmaDomain::pos_x)(i,0));
  for(int j=0; j<m_pd.m_ydim; j++) y_vec.push_back(m_pd.grid(PlasmaDomain::pos_y)(0,j));
  if(m_pd.m_continue_mode){ // Check for existing end.tpstate file when in continue mode
    m_init_filename = m_pd.m_out_directory/fs::path("end.tpstate");
    fs::directory_entry init_path_dir(m_init_filename);
    assert(init_path_dir.exists() && init_path_dir.is_regular_file() && "Continue directory must contain end.tpstate file for tracer particles");
  } else { // Copy init.tpstate into output directory, if not continue mode
    // read the init tp state file here
    fs::path new_init_path = m_pd.m_out_directory/"init.tpstate";
    fs::directory_entry new_init_dir(new_init_path);
    if(new_init_dir.exists()) fs::remove(new_init_path);
    fs::copy(m_init_filename, new_init_path, fs::copy_options::overwrite_existing);
    fs::path old_out_path = m_pd.m_out_directory/"particles.tpout";
    fs::directory_entry old_out_dir(old_out_path);
    if(old_out_dir.exists()) fs::remove(old_out_path);
  }
  //Either way, read in init tpstate file
  readTPStateFile(m_init_filename);
  if(!m_pd.m_continue_mode){
    writeToTPOutFile(m_pd.m_out_directory/m_out_filename, 0.0);
  }
}

void TracerParticles::parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs){
  for(int i=0; i<lhs.size(); i++){
    std::string this_lhs = lhs[i];
    std::string this_rhs = rhs[i];
    if(this_lhs == "init_file") m_init_filename = fs::path(this_rhs);
    // else if(this_lhs == "duration") duration = std::stod(this_rhs);
    else std::cerr << this_lhs << " config not recognized.\n";
  }
}

void TracerParticles::iterateModule(double dt){
  Grid &m_v_x = m_pd.m_eqs->grid(IdealMHD::v_x);
  Grid &m_v_y = m_pd.m_eqs->grid(IdealMHD::v_y);
  for(int i=m_particles.size()-1; i>=0; i--){
    std::vector<double> &p = m_particles[i];
    double v_x_p = bilinearInterpolate(p, m_v_x, x_vec, y_vec);
    double v_y_p = bilinearInterpolate(p, m_v_y, x_vec, y_vec);
    std::vector<double> p_half_step = {p[0] + 0.5*dt*v_x_p, p[1] + 0.5*dt*v_y_p};
    double v_x_p_half = bilinearInterpolate(p_half_step, m_v_x, x_vec, y_vec);
    double v_y_p_half = bilinearInterpolate(p_half_step, m_v_y, x_vec, y_vec);
    p[0] += dt*v_x_p_half; p[1] += dt*v_y_p_half;
    if(p[0] < x_vec[0] || p[0] > x_vec.back()){
      if(m_pd.x_bound_1 == PlasmaDomain::BoundaryCondition::Periodic){
        double width = (x_vec.back() - x_vec[0]);
        p[0] = std::fmod((p[0] - x_vec[0] + width),width) + x_vec[0];
      } else {
        m_particles.erase(m_particles.begin() + i);
      }
    }
    if(p[1] < y_vec[0] || p[1] > y_vec.back()){
      if(m_pd.y_bound_1 == PlasmaDomain::BoundaryCondition::Periodic){
        double height = (y_vec.back() - y_vec[0]);
        p[1] = std::fmod((p[1] - y_vec[0] + height),height) + y_vec[0];
      } else {
        m_particles.erase(m_particles.begin() + i);
      }
    }
  }
  int old_time_iter = (int)(m_pd.m_time/m_pd.m_time_output_interval);
  int new_time_iter = (int)((m_pd.m_time + dt)/m_pd.m_time_output_interval);
  bool store_cond_1 = m_pd.m_iter_output_interval > 0 && (m_pd.m_iter+1)%m_pd.m_iter_output_interval == 0;
  bool store_cond_2 = m_pd.m_time_output_interval > 0.0 && new_time_iter > old_time_iter;
  if (store_cond_1 || store_cond_2) {
    writeToTPOutFile(m_pd.m_out_directory/m_out_filename, dt);
  }
  writeTPStateFile(m_pd.m_out_directory/m_end_filename);
}


std::string TracerParticles::commandLineMessage() const
{
  return "Tracer Particles On";
}


void TracerParticles::readTPStateFile(const fs::path &init_path) {
  fs::directory_entry init_path_dir(init_path);
  assert(init_path_dir.exists() && init_path_dir.is_regular_file() && "Tracer particle initialization file must exist");
  std::ifstream init_file(init_path.string());
  std::string line;
  // read through config file
  while (std::getline(init_file, line)){
    // obtain the lhs and rhs quantities for the current line
    clearWhitespace(line);
    if (line.empty() || line[0] == '#') continue; // skip comment and empty lines
    std::cout << line << std::endl;
    std::vector<std::string> components_str = splitString(splitString(line,'#')[0],',');
    m_particles.push_back({std::stod(components_str[0]),std::stod(components_str[1])});
  }
  init_file.close();
}

void TracerParticles::writeTPStateFile(const fs::path &state_path) {
  std::ofstream state_file((m_pd.m_out_directory/m_end_filename).string());
  state_file.precision(std::numeric_limits<double>::digits10 + 1);
  for(std::vector<double> p : m_particles){
    state_file << p[0] << "," << p[1] << std::endl;
  }
  state_file.close();
}

void TracerParticles::writeToTPOutFile(const fs::path &out_path, double dt) {
  std::ofstream out_file((m_pd.m_out_directory/m_out_filename).string(),std::ofstream::app);
  out_file.precision(std::numeric_limits<double>::digits10 + 1);
  out_file << "t=" << m_pd.m_time + dt << std::endl;
  for(std::vector<double> p : m_particles){
    out_file << p[0] << "," << p[1] << std::endl;
  }
  out_file.close();
}

// std::vector<std::vector<double>> TracerParticles::readParticleCoordinates(std::string s){
//   std::string particle;
//   std::istringstream s_particles(s);
//   while(getCleanedLine(s_particles, particle, ';')){
//     double x, y;
//     std::string coord;
//     std::istringstream s_particle(particle);
//     std::getline(s_particle,coord,',');
//     double x = std::stod(coord);
//     std::getline(s_particle,coord,',');
//     double y = std::stod(coord);
//     m_particles.push_back({x,y});
//   }
// }
