function [lambda] = getCoulombLogarithm(n,T)

Gam_e = getGamma(n,T);
lambda = 1./sqrt(3.*Gam_e.^3);

end