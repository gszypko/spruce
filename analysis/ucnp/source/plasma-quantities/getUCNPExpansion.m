function [sig,Te,v] = getUCNPExpansion(t,sig0,Te0,x,y)
% t (vector double): time in s
% sig0 (double): initial geometric mean of rms plasma size in cm
% Te0 (double): initial electron temperature in K
% x (vector double): x position in cm
% y (vector double): y position in cm

tau = getTauExp(sig0,2*Te0);
sig = sig0.*sqrt(1+t.^2./tau^2);
Te = Te0./(1+t.^2./tau^2);

v = struct;
for i = 1:length(t)
    gam = t(i)./tau^2./(1+t(i).^2/tau^2);
    v(i).vx = gam.*x;
    v(i).vy = gam.*y;
end


end