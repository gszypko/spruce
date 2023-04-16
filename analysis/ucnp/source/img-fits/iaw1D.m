function [f] = iaw1D(x,iaw_amp,iaw_lamda,iaw_sig)
G = @(x,sig) exp(-(x).^2./(2*sig.^2));
P = @(x,sig) -iaw_amp*cos(2*pi.*x./iaw_lamda).*G(x,sig);
% P2 = @(x,sig) iaw_amp_2*cos(2*2*pi.*x./iaw_lamda).*G(x,sig);
f = P(x,iaw_sig);
end