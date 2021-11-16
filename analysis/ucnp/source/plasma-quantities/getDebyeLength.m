function [lambda] = getDebyeLength(n,Te)
% n (mat double): plasma number density (cm^-3)
% Te (mat double): electron temperature (K)
% lambda (mat double): electron debye screening length (cm)

c = defineConstants(); 
lambda = sqrt(c.kB.*Te./(4.*pi.*n.*c.e^2));

end