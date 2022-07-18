clearvars
setpath;
c = defineConstants();
n = 1e9;
Te = 20;
v_th = sqrt(c.kB*Te/c.mE);
total_width = .6;
cell_width = total_width/101;
tau = cell_width/(3*v_th);
w_pe = getPlasmaFreq(n,c.mE);
tau_pe = 1/w_pe;