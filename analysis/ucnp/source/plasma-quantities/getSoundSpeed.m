function [cs] = getSoundSpeed(m_i,Ti,Te,gam_e,gam_i)
% m_i (double): ion mass
% Te (double): electron temperature
% gam (double): adiabatic index
% Note: all units cgs
cs = sqrt(cts.cgs.kB*(gam_i*Ti+gam_e*Te)/m_i);
end