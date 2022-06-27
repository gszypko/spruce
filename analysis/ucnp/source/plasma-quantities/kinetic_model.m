function [t,sig,gam,Ti,Te,n0] = kinetic_model(t0,sig0,Ti0,Te0,n0)
% solves the kinetic equations from https://aip.scitation.org/doi/10.1063/1.4915135
% includes electron-ion thermalization but not disorder-induced heating term

% define constants and function for thermalization rate
c = defineConstants();
kB = c.kB; m_i = c.mI; m_e = c.mE;
nu_ei = @(n,Te) 2*(m_e/m_i)*getEICRate(n,Te,m_e);
if nu_ei(n0,Te0) < 0, error("Electrons are strongly coupled."); end

% define rate equations for ode45 solver
    function [dydt] = rate_equations(~,y)
        % y = [sig^2; gam; Ti; Te]
        n = n0*sig0^3/y(1)^(3/2);
        dydt = [0; 0; 0; 0];
        % square plasma size
        dydt(1) = 2*y(2)*y(1);
        % hydrodynamic expansion parameter
        dydt(2) = kB*(y(3)+y(4))/(m_i*y(1)) - y(2)^2;
        % ion temperature
        dydt(3) = -2*y(2)*y(3) + nu_ei(n,y(4))*(y(4)-y(3));
        % electron temperature
        dydt(4) = -2*y(2)*y(4) - nu_ei(n,y(4))*(y(4)-y(3));
    end

% do ode solving
y0 = [sig0^2; 0; Ti0; Te0];
options = odeset('RelTol',1e-5);
[t,y] = ode45(@(t,y) rate_equations(t,y),t0,y0,options);

% output results
sig = sqrt(y(:,1));
gam = y(:,2);
Ti = y(:,3);
Te = y(:,4);
n = n0./sig.^3;

end