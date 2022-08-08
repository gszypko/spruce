function [lambda] = getDebyeLength(n,Te)
% n (mat double): plasma number density (cm^-3)
% Te (mat double): electron temperature (K)
% lambda (mat double): electron debye screening length (cm)

lambda = sqrt(cts.cgs.kB.*Te./(4.*pi.*n.*cts.cgs.e^2));

end