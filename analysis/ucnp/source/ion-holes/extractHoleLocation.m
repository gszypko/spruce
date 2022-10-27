function [fit] = extractHoleLocation(t,r,n,v,sig0,tau,r_hole)
% t (double): current time
% r (vector): position along hole propagation axis
% n (vector): ion density for each position 'r'
% v (vector): fluid velocity along propagation axis for each position 'r'
% sig0 (double): initial geometric mean of rms plasma size
% tau (double): hydrodynamic expansion timescale
% r_hole (double): expected hole position relative to plasma center for given time point

% compute various things
sig = sig0*sqrt(1+t^2/tau^2); % expected geometric mean of rms plasma size
gam = t/tau^2/(1+t^2/tau^2); % hydrodynamic expansion parameter (i.e., slope of velocity profile)
pi = n.*v; % momentum density (in this case does not have the mass...)

% filter out any potential noise from n and v
n_sg = sgolayfilt(n,3,5);
v_sg = sgolayfilt(v,3,5);

% define general anonymous functions
G = @(r,A,r0,sig) A.*exp(-(r-r0).^2./(2*sig^2));
L = @(r,m,b) m.*r + b;

% fit gaussian to density profile
fun = @(p,r) G(r,p(1),p(2),p(3));
p0(1) = max(n_sg); lb(1) = min(n_sg); ub(1) = max(n_sg)*2;
p0(2) = 0; lb(2) = min(r); ub(2) = max(r);
p0(3) = sig; lb(3) = sig/5; ub(3) = sig*5;
p = lsqcurvefit(fun,p0,r,n,lb,ub);
g_fit.n_max = p(1);
g_fit.r0 = p(2);
g_fit.sig = p(3);

% fit linear model to velocity profile
fun = @(p,r) L(r,p(1),p(2));
p0(1) = gam; lb(1) = gam/5; ub(1) = gam*5;
p0(2) = 0; lb(2) = -max(abs(v)); ub(2) = max(abs(v));
p = lsqcurvefit(fun,p0,r,v,lb,ub);
linear_fit.gam = p(1);
linear_fit.v0 = p(2);

% normalize quantities for constrained fit
n_norm = mean(abs(n));
v_norm = mean(abs(v));
if v_norm == 0, v_norm = 1; end
n_tild = n./n_norm;
v_tild = v./v_norm;

% fit model for hole Location
    function [model] = constrained_model(p,r)
        % unfold position vector: r = [r_g r_v r_pi]
        
        % create fit models with perturbations
        Gp = @(r,A,r0,sig,Ag,sigp,rL,rR,Ac) G(r,A,r0,sig) - G(r,Ag,rL,sigp) - G(r,Ag,rR,sigp) - G(r,Ac,0,sigp);

        % create full concotenated fit model
            % p = [A r0 sig Ag sigp rL rR m b Av Ac rL2 rR2]
        model = Gp(r,p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8));
    end

% constrained fit routine
p0(1) = max(n)/n_norm*1.05; lb(1) = p0(1)/1.25; ub(1) =p0(1)*1.25;
p0(2) = 0; lb(2) = -sig0/4; ub(2) =sig0/4;
p0(3) = g_fit.sig; lb(3) = sig/1.25; ub(3) = sig*1.25;
p0(4) = .05*max(g_fit.n_max)/n_norm; lb(4) = 0; ub(4) = p0(4)*10;
p0(5) = .15*max(g_fit.sig); lb(5) = p0(5)/10; ub(5) = p0(5)*10;
p0(6) = -r_hole; lb(6) = min(r); ub(6) = 0;
p0(7) = r_hole; lb(7) = 0; ub(7) = max(r);
p0(8) = .1*max(g_fit.n_max)/n_norm; lb(8) = 0; ub(8) = p0(8)*10;

fun = @(p,r) constrained_model(p,r);
xdata = [r];
ydata = [n_tild];
p = lsqcurvefit(fun,p0,xdata,ydata,lb,ub);

% output fit parameters
fit.t = t;
fit.g_amp = p(1)*n_norm;
fit.g_r0 = p(2);
fit.g_sig = p(3);
fit.p_xL = p(6);
fit.p_xR = p(7);
fit.p_n_amp = p(4)*n_norm;
fit.p_sig = p(5);


% fig = figure;
% plot(xdata,ydata,'.')
% hold on
% plot(xdata,fun(p0,xdata),'.')
% plot(xdata,fun(p,xdata),'.')
% legend({'data','guess','fit'})
% pause(1)
% close(fig)
end