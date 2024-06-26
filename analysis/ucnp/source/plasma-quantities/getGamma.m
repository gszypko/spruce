function [Gam] = getGamma(n,T)
% n (mat double): plasma number density (cm^-3)
% T (mat double): temperature of plasma species (K) - ion or electron
% Gam (mat double): Coulomb coupling parameter (dimensionless)

a = getWignerSeitzRadius(n); % average interparticle spacing (cm)
Gam = cts.cgs.e^2./(cts.cgs.kB.*T.*a); % calculate Coulomb coupling parameter

end