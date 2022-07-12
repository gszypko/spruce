function [beta] = getPlasmaBeta(n,T,B)
% n (nxm double): plasma number density (cm^-3)
% T (nxm double): plasma temperature (K)
% B (nxm double): magnetic field amplitude (G)
% beta (nxm double): ratio of plasma pressure (nkT) to magnetic pressure (B^2/8pi)
% Note - all units in cgs
c = defineConstants();
beta = n.*c.kB.*T./(B.^2./(8*pi));

end