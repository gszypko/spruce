function [r_c] = getGyroRadius(B,T,m)
% all quantities are in cgs units
% B (nxm double): magnetic field amplitude
% T (nxm double): temperature of given species
% m (double): mass of given species
% r_c: cyclotron (gyro) radius
c = defineConstants();
v_th = sqrt(c.kB.*T./m);
w_c = getGyroFreq(B,m);
r_c = v_th./w_c;
end