function [tau] = getTauExp(sig,Te0,usebeta,beta)
% sig (mat double): geometric mean of RMS plasma size (cm)
% Te0 (mat double): initial electron temperature (K)
% usebeta (bool): (true) use cuspy defintion (false) use gaussian definition
% beta (double): mo
% tau_exp (mat double): hydrodynamic expansion timescale (s)

c = defineConstants();
if nargin < 3, usebeta = false; beta = 1; end
if usebeta && nargin < 4, beta = 0.63; end
tau = sqrt(c.mI.*sig.^2./(c.kB.*Te0)).*beta;

end