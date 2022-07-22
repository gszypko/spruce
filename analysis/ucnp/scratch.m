clearvars
setpath;
c = defineConstants();

sig = 0.1;
total_width = 10*sig;
cell_width = total_width/101;

n = 1e9;
Te = 20;
B = 10;

l_deb = getDebyeLength(n,Te);
k_small = 2*pi/(cell_width);
k_large = 2*pi/total_width;


w_pe = getPlasmaFreq(n,c.mE);
tau_pe = 1/w_pe;

w_c = getGyroFreq(B,c.mE);

v_th = sqrt(c.kB*Te/c.mE);
w_th_small = sqrt(3)*v_th*k_small;
w_th_large = sqrt(3)*v_th*k_large;

w_tot_small = sqrt(w_pe^2+w_th_small^2);
w_tot_large = sqrt(w_pe^2+w_th_large^2);

v_pe = w_tot_small/k_small;

dt = cell_width/(v_pe);