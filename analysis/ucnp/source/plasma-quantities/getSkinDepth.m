function [d_s] = getSkinDepth(n)
% n (nxm double): plasma number density (cm^-3)
% d_s (nxm double): plasma skin depth (cm)
    % this is the distance that EM radiation can penetrate into a plasma
w_pe = getPlasmaFreq(n,c.mE);
d_s = cts.cgs.c./w_pe;
end