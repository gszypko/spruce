function [w_c] = getGyroFreq(B,m)
% B (nxm double): magnetic field amplitude (G)
% m (double): particle mass (g)
% w_c (nxm double): gyro frequency for particle of mass m (s^-1)
w_c = cts.cgs.e*B/(m*cts.cgs.c);
end