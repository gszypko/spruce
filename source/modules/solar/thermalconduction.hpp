//thermalconduction.hpp
//Header for the Thermal Conduction Module,
//an implementation of the abstract Module class
//Applies field-aligned Spitzer-Harm thermal conductivity
//with free-streaming saturation

#ifndef THERMALCONDUCTION_HPP
#define THERMALCONDUCTION_HPP

#include "module.hpp"

class PlasmaDomain;

class ThermalConduction : public Module {
    public:
        ThermalConduction(PlasmaDomain &pd);
        void parseModuleConfigs(std::vector<std::string> lhs, std::vector<std::string> rhs);
        void preIterateModule(double dt) override;
        void iterateModule(double dt) override;
        void setupModule() override;
        void fileOutput(std::vector<std::string>& var_names, std::vector<Grid>& var_grids) override;
        std::string commandLineMessage() const override;
    private:
        bool flux_saturation;
        double epsilon;
        double dt_subcycle_min;
        int curr_num_subcycles;
        bool output_to_file = false;
        Grid avg_change;
        std::string time_integrator;
        bool inactive_mode = false; //when true, module will compute the change in thermal energy (allowing for output), but WON'T apply it to the simulation
        int numberSubcycles(double dt);
        Grid oneDimConductiveFlux(const Grid &temp, const Grid &rho, double k0, int index) const;
        void fieldAlignedConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &temp, const Grid &rho,
                                            const Grid &b_hat_x, const Grid &b_hat_y, double k0) const;
        void saturateConductiveFlux(Grid &flux_out_x, Grid &flux_out_y, const Grid &rho, const Grid &temp) const;
        Grid saturatedKappa() const;
        std::vector<Grid> saturationTerms() const;
        Grid thermalEnergyDerivative(Grid &m_temp, std::vector<Grid> b_hat) const;
};

#endif