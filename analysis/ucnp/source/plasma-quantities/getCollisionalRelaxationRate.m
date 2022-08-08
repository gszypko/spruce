function [mu] = getCollisionalRelaxationRate(n,Ti,a,b)
% n (double): plasma number density (cm^-3)
% Ti (double): ion temperature (K)
% a (double): parameterizes collisional relaxation rate function (dimensionless)
% b (double): parameterizes collisional relaxation rate function (dimensionless)

if nargin < 3 % default values for a and b if not specified
    a = 1;
    b = 0.7;
elseif nargin == 3 % default value for b if not specified
    b = 0.7;
end

Gam = getGamma(n,Ti); % ion coulomb coupling parameter
wp = getPlasmaFreq(n,cts.cgs.mI); % ion plasma frequency with units gam422
Lam = 1/sqrt(3*Gam^3); % plasma parameter
mu = a*Gam^(3/2)*wp*log(1+b*Lam); % collisional relaxation rate with units gam422 (see PhysRevLett.110.235001)

end

