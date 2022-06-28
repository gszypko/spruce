function [gam_ei] = getEICRate(n,Te,m_e)
Gam_e = getGamma(n,Te);
w_pe = getPlasmaFreq(n,m_e);
LamE = getPlasmaParameter(n,Te);
gam_ei = sqrt(2/3/pi).*Gam_e.^(3/2).*w_pe.*log(LamE);

end