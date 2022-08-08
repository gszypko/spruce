function [r_L] = getLandauLength(T)
% Landau length is the classical distance of approach
% All units cgs
% T is the temperature of the plasma species in question
r_L = cts.cgs.e^2./(cts.cgs.kB.*T);
end