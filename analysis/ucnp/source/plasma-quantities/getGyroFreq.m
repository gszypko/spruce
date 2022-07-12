function [w_c] = getGyroFreq(B,m)
% B (nxm double): magnetic field amplitude (G)
% m (double): particle mass (g)
% w_c (nxm double): gyro frequency for particle of mass m (s^-1)
c = defineConstants();
w_c = c.e*B/(m*c.c);
end