function [results] = fitGaussianNoOffset1D(x,y)
% x (1xn double): independent data
% y (1xn double): dependent data corresponding to x, same length as x

% apply savitzky-golay filter to y data for determining initial guesses
polynomial_order = 3;
kernel_size = 11;
y_sg = sgolayfilt(y,polynomial_order,kernel_size);

% initial guesses for fit
amp = max(y_sg);
cen = trapz(x,x.*y_sg)/trapz(x,y_sg);
sig = sqrt(trapz(x,(x-cen).^2.*y_sg)/trapz(x,y_sg));

% define fit model
G = @(x,A,x0,sig) A.*exp(-(x-x0).^2./(2*sig^2));
fit_model = @(p,x) G(x,p(1),p(2),p(3));

% initialize parameter bounds
p0 = []; lb = []; ub = [];
p0(1) = amp; lb(1) = 0; ub(1) = amp*10;
p0(2) = cen; lb(2) = min(x); ub(2) = max(x);
p0(3) = sig; lb(3) = 0; ub(3) = sig*10;

% do fit and evaluate error bounds
p = lsqcurvefit(fit_model,p0,x,y,lb,ub);
% pci = nlparci(p,R,'jacobian',J);    % 95% confidence intervals from fit
% pse = (pci(:,2) - pci(:,1))/3.92;    % convert 95% confidence intervals to standard error

% output results
results.fit = fit_model(p,x);
results.guess = fit_model(p0,x);
results.amp = p(1);
results.cen = p(2);
results.sig = p(3);

end