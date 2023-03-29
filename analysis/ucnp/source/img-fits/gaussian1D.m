function [f] = gaussian1D(x,amp,cen,sig)
f = amp.*exp(-(x-cen).^2./(2*sig.^2));
end