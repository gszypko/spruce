function [cs] = getSoundSpeed(m_i,Te,gam)
% m_i (double): ion mass
% Te (double): electron temperature
% gam (double): adiabatic index
% Note: all units cgs
cs = sqrt(gam*cts.cgs.kB*Te/m_i);
end