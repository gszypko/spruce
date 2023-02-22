function [tau] = getTauExp(sig,T0,usebeta,beta)
% sig (mat double): geometric mean of RMS plasma size (cm)
% T0 (mat double): initial plasma temperature (K)
% usebeta (bool): (true) use cuspy defintion (false) use gaussian definition
% beta (double): mo
% tau_exp (mat double): hydrodynamic expansion timescale (s)

if nargin < 3, usebeta = false; beta = 1; end
if usebeta && nargin < 4, beta = 0.63; end
tau = sqrt(cts.cgs.mI.*sig.^2./(cts.cgs.kB.*T0)).*beta;

end