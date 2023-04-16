function [results] = fitIAWDispersion(t,A,w0,tau)

    function [amp] = amp_fun(t,A0,g0,w0,tau)
        temp = @(t) w0.*integral(@(tp) 1./(1+tp.^2./tau^2),0,t);
        phi = zeros(size(t));
        for i1 = 1:length(t)
            phi(i1) = temp(t(i1));
        end
        amp = A0.*exp(-g0.*t).*cos(phi);
    end
fit_model = @(c,d) amp_fun(d,c(1),c(2),c(3),c(4));
p0 = []; lb = []; ub = [];
p0(1) = A(1); lb(1) = -10*A(1); ub(1) = 10*A(1);
p0(2) = 0.2/tau; lb(2) = 0; ub(2) = 10*p0(2);
p0(3) = w0; lb(3) = w0/10; ub(3) = 10*w0;
p0(4) = tau; lb(4) = tau; ub(4) = tau;

% do fitting and error estimation
p = lsqcurvefit(fit_model,p0,t,A,lb,ub);
% pci = nlparci(p,R,'jacobian',J);    % 95% confidence intervals from fit
% pse = (pci(:,2) - pci(:,1))/3.92;    % convert 95% confidence intervals to standard error

% output results
results.fit = fit_model(p,t);
results.guess = fit_model(p0,t);
results.A0 = p(1);
results.gam = p(2);
results.w0 = p(3);
results.tau = p(4);

end