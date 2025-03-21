function [p,y_fit,y_guess,pci,pse,fit_model] = fitSlopedGaussianWithIAW1D(x,y,iaw_amp,iaw_lam)

y_sg = sgolayfilt(y,3,9);
amp = max(y_sg);
x_cen = trapz(x,x.*y_sg)/trapz(x,y_sg);
sig = sqrt(trapz(x,(x-x_cen).^2.*y_sg)/trapz(x,y_sg));

fit_model = @(p,x) slopedGaussianWithIAW1D(x,p(1),p(2),p(3),p(4),p(5),p(6),p(7),p(8),p(9),p(10));
p0 = []; lb = []; ub = [];
p0(1) = amp*0.9; lb(1) = 0; ub(1) = amp*10;
p0(2) = x_cen; lb(2) = min(x); ub(2) = max(x);
p0(3) = sig; lb(3) = sig./10; ub(3) = sig*10;
p0(4) = iaw_amp; lb(4) = .01; ub(4) = 0.5;
p0(5) = iaw_lam; lb(5) = p0(5)/10; ub(5) = p0(5)*10;
p0(6) = 0; lb(6) = -pi; ub(6) = pi;
p0(7) = sig; lb(7) = p0(7)/10; ub(7) = p0(7)*10;
p0(8) = 0; lb(8) = -inf; ub(8) = inf;
p0(9) = 0; lb(9) = min(x); ub(9) = max(x);
p0(10) = 0.01; lb(10) = .001; ub(10) = 0.1;

[p,~,R,~,~,~,J] = lsqcurvefit(fit_model,p0,x,y,lb,ub);

pci = nlparci(p,R,'jacobian',J);    % 95% confidence intervals from fit
pse = (pci(:,2) - pci(:,1))/3.92;    % convert 95% confidence intervals to standard error

y_fit = fit_model(p,x);
y_guess = fit_model(p0,x);

end