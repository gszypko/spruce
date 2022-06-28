n0 = 1e9;
Ti0 = 1;
Te0 = 20;
sig0 = .1;
tau = getTauExp(sig0,Te0+Ti0);
t = linspace(0,2*tau,101);
[t,sig,gam,Ti,Te,n0] = kinetic_model(t,sig0,Ti0,Te0,n0);

figure
plot(t./tau,Ti)