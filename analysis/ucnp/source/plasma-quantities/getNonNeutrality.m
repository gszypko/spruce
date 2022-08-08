function [frac] = getNonNeutrality(n0,Te,sig)
% n0 (mat double): peak plasma density for exponential plasma (cm^-3)
% Te (mat double): initial electron temperature (K)
% sig (mat double): geometric mean of rms plasma size (cm)

% See page 17 of Thomas Langin's PhD thesis:
% https://ultracold.rice.edu/publications/ThomasLanginPhDThesis.pdf

spacelim = sig*100;
x = linspace(-spacelim,spacelim,5001);
y = linspace(-spacelim,spacelim,5001);
z = linspace(-spacelim,spacelim,5001);
n = @(x,y,z) n0.*exp(-sqrt(x.^2+y.^2/4+z.^2/4)./sig);

Ni = integral3(n,min(x),max(x),min(y),max(y),min(z),max(z)); % number of ions
Ns = 3/2*sqrt(pi/2)*(4*pi*cts.si.eps*sig/cts.si.e^2)*c.kB*Te; % N^star from Thomas's thesis
frac = (sqrt(Ni/Ns)-1)/sqrt(Ni/Ns); % fraction of trapped electrons

end