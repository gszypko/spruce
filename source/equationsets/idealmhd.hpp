//idealmhd.hpp
//Defines IdealMHD EquationSet

#ifndef IDEALMHD_HPP
#define IDEALMHD_HPP

#include "grid.hpp"
#include "equationset.hpp"
#include <vector>
#include <string>

class PlasmaDomain;

class IdealMHD: public EquationSet {
    public:
        IdealMHD(PlasmaDomain &pd);
        IdealMHD() = default;
        
        std::vector<std::string> def_var_names() const override{
            return {"rho","temp","mom_x","mom_y","mom_z","bi_x","bi_y","bi_z","grav_x","grav_y",
                "n","press","thermal_energy","v_x","v_y","v_z","kinetic_energy",
                "b_x","b_y","b_z","b_mag","b_hat_x","b_hat_y","b_hat_z","dt"};
        }
        enum Vars {rho,temp,mom_x,mom_y,mom_z,bi_x,bi_y,bi_z,grav_x,grav_y,
                n,press,thermal_energy,v_x,v_y,v_z,kinetic_energy,
                b_x,b_y,b_z,b_mag,b_hat_x,b_hat_y,b_hat_z,dt};

        std::vector<int> state_variables() const override {
            return {rho,temp,mom_x,mom_y,mom_z,bi_x,bi_y,bi_z,grav_x,grav_y};
        }

        std::vector<int> evolved_variables() const override {
            return {rho,mom_x,mom_y,mom_z,thermal_energy,bi_x,bi_y,bi_z};
        }

        std::vector<int> densities() const override { return {rho}; }
        std::vector<std::vector<int>> momenta() const override { return {{mom_x,mom_y}}; }
        std::vector<int> thermal_energies() const override { return {thermal_energy}; }
        std::vector<int> fields() const override { return {bi_x,bi_y}; }
        std::vector<int> timescale() const override {return {dt}; }

    private:
        double m_global_viscosity{0};
        std::string m_viscosity_opt{"velocity"};
        std::vector<Grid> computeTimeDerivativesDerived(const std::vector<Grid>& grids) const override;
        void setupEquationSetDerived() override;
        void recomputeEvolvedVarsFromStateVars(std::vector<Grid> &grids) const override;
        void recomputeDerivedVarsFromEvolvedVars(std::vector<Grid> &grids) const override;
        void catchNullFieldDirection(std::vector<Grid> &grids) const;
        void recomputeDT(std::vector<Grid>& grids) const override;
        std::vector<Grid> computeTimeDerivativesCharacteristicBoundary(const std::vector<Grid> &grids, bool x_bound_1, bool x_bound_2, bool y_bound_1, bool y_bound_2, double visc_coeff) const;
        std::vector<std::vector<Grid>> singleBoundaryTermsMOC(const std::vector<Grid> &grids, int boundary_index, bool boundary_lower, double visc_coeff) const;
        void parseEquationSetConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs) override;
        bool moc_b_limiting{false};
        double moc_b_lower_lim{0.1}; //Smallest fraction of nearest interior value allowed in MoC ghost zone. Larger values lead to stricter thresholding. Should be <=1.0, can be negative.
        double moc_b_upper_lim{10.0}; //Largest multiple of nearest interior value allowed in MoC ghost zone. Larger values lead to stricter thresholding. Should be >=1.0.
        bool moc_mom_limiting{false};
        double moc_mom_lower_lim{0.1}; //Smallest fraction of nearest interior value allowed in MoC ghost zone. Larger values lead to stricter thresholding. Should be <=1.0, can be negative.
        double moc_mom_upper_lim{10.0}; //Largest multiple of nearest interior value allowed in MoC ghost zone. Larger values lead to stricter thresholding. Should be >=1.0.
        void applyBThresholdingMoC(std::vector<Grid> &grids) const;
        void applyMomThresholdingMoC(std::vector<Grid> &grids) const;
};

#endif