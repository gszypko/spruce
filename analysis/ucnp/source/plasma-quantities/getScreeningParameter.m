function [kappa] = getScreeningParameter(n,Te)
% n (mat double): plasma number density (cm^-3)
% Te (mat double): electron temperature (K)
% kappa (mat double): electron screening parameter (dimensionless)

a = getWignerSeitzRadius(n);
lambda = getDebyeLength(n,Te);
kappa = a./lambda;

end