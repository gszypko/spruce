function [mu] = getCollisionalRelaxationRateMD(n,Ti)
% n (double): ion density (cm^-3)
% Ti (double): ion temperature (K)

% get Coulomb coupling parameter and plasma frequency corresponding to user input
Gam = getGamma(n,Ti); % ion coulomb coupling parameter
wp = getPlasmaFreq(n,cts.cgs.mI); % ion plasma frequency with units gam422

% tabulated Gamma and mu values from our 2012 velocity relaxation paper (https://link.aps.org/doi/10.1103/PhysRevLett.109.185008)
% these were copy and pasted from '\\10.65.11.22\Killian_Research\Plasma\vccAnalysis\digitized_paper_data', which is on the
% drobo. This is the MD data that is plotted in Fig. 5 of 2012 paper (i.e., thick solid blue line with dots), which gives the
% collisional relaxation rate per unit plasma frequency that is correct even in the strongly coupled regime.
Gam_tab = [1.3215e-02 2.1451e-02 4.6388e-02 9.9675e-02 2.1692e-01 4.6611e-01 1.0080e+00 2.1521e+00 4.6539e+00 1.0064e+01];
mu_tab = [3.6831e-03 6.6173e-03 1.5581e-02 3.3525e-02 6.0234e-02 1.0191e-01 1.6359e-01 2.3638e-01 3.2164e-01 4.3113e-01];

% this function interpolates the tabulated data to return the collisional relaxation rate that corresponds to the user
% supplied density and temperature values.
mu = wp*interp1(Gam_tab,mu_tab,Gam,'makima','extrap');
end