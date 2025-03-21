function [f] = gaussianWithIAW1D(x,amp,cen,sig,iaw_amp,iaw_lamda,iaw_phase,iaw_sig,iaw_amp_2)
G = @(x,sig) amp.*exp(-(x-cen).^2./(2*sig.^2));
P = @(x,sig) -iaw_amp*cos(2*pi.*x./iaw_lamda + iaw_phase).*G(x,sig);
P2 = @(x,sig) iaw_amp_2*cos(2*2*pi.*x./iaw_lamda + iaw_phase).*G(x,sig);
f = G(x,sig) + P(x,iaw_sig);
end