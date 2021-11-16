function [w_p] = getPlasmaFreq(n,m)
% n (mat double): plasma number density (cm^-3)
% m (mat double): particle mass (kg)
% w_p (mat double): plasma oscillation frequency (rad/s)

c = defineConstants();
w_p = sqrt(4.*pi.*n.*c.e^2./m);

end