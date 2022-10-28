function [r] = ionHoleProp(t,v,Ti,Te,m_i,gam)
% t (vector): time points to return hole position for
% sig (double): initial geometric mean of plasma size
% tau (double): hydrodynamic expansion timescale
% v (fun): v(t,r) returns velocity at position r and time t
% Ti (fun): Ti(t,r) returns velocity at position r and time t
% Te (fun): Te(t,r) returns velocity at position r and time t
% m_i (double): ion mass
% gam (double): adiabatic index

% speed of sound Cs = dP/drho (i.e., change in total pressure w.r.t. mass density)
% 'gam' can be varied depending on the adiabaticity of the plasma

Cs = @(t,r) sqrt(gam*cts.cgs.kB*(Te(t,r)+Ti(t,r))/m_i);
odefun = @(t,r) v(t,r) + Cs(t,r);
r0 = 0;
if length(t) == 1
    tspan = [0 t];
else
    tspan = t;
end
if length(t) > 1 || t(end) > 0
    [~,rsol] = ode45(odefun,tspan,r0);
else
    rsol = 0;
end
if length(t)==1
    r = rsol(end);
else
    r = rsol;
end

end