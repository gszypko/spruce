function [p,y_fit,y_guess,pci,pse,fit_model] = fitGaussian1D(x,y)

y_sg = sgolayfilt(y,3,9);
amp = max(y_sg);
x_cen = trapz(x,x.*y_sg)/trapz(x,y_sg);
sig = sqrt(trapz(x,(x-x_cen).^2.*y_sg)/trapz(x,y_sg));

fit_model = @(p,x) gaussian1D(x,p(1),p(2),p(3));
p0 = []; lb = []; ub = [];
p0(1) = amp; lb(1) = 0; ub(1) = amp*10;
p0(2) = x_cen; lb(2) = min(x); ub(2) = max(x);
p0(3) = sig; lb(3) = sig./10; ub(3) = sig*10;

[p,~,R,~,~,~,J] = lsqcurvefit(fit_model,p0,x,y,lb,ub);

pci = nlparci(p,R,'jacobian',J);    % 95% confidence intervals from fit
pse = (pci(:,2) - pci(:,1))/3.92;    % convert 95% confidence intervals to standard error

y_fit = fit_model(p,x);
y_guess = fit_model(p0,x);

end