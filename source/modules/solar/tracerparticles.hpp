//tracerparticles.hpp
//Header for the Tracer Particles Module,
//an implementation of the abstract Module class
//Tracks a set of Lagrangian point particles that follow plasma flow

#ifndef TRACERPARTICLES_HPP
#define TRACERPARTICLES_HPP

#include "module.hpp"
#include <vector>
#include <filesystem>
namespace fs = std::filesystem;

class PlasmaDomain;

class TracerParticles : public Module {
    public:
        TracerParticles(PlasmaDomain &pd);
        void iterateModule(double dt) override;
        std::string commandLineMessage() const override;
        void setupModule() override;
    private:
        fs::path m_out_filename{"particles.tpout"}; // file name for saving time evolution of tracer particles
        fs::path m_init_filename{"init.tpstate"};
        fs::path m_end_filename{"end.tpstate"};
        std::vector<double> x_vec, y_vec;
        std::vector<std::vector<double>> m_particles;
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        void readTPStateFile(const fs::path &init_path);
        void writeTPStateFile(const fs::path &state_path);
        void writeToTPOutFile(const fs::path &out_path, double dt);
};

#endif