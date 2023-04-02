function [f] = slopedGaussian1D(x,amp,cen,sig,slope,slope_center)
f = amp.*exp(-(x-cen).^2./(2*sig.^2)) + slope.*(x-slope_center);
end