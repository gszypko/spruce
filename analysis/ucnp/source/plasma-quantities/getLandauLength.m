function [r_L] = getLandauLength(T)
% Landau length is the classical distance of approach
% All units cgs
% T is the temperature of the plasma species in question
c = defineConstants();
r_L = c.e^2./(c.kB.*T);
end