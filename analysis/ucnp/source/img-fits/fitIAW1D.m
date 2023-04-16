function results = fitIAW1D(x,y,k,offset)
% x (1xn double): independent data
% y (1xn double): dependent data corresponding to x, same length as x

if nargin >= 4, offset = mean(y); 
else offset = 0; end


% apply savitzky-golay filter to y data for determining initial guesses
polynomial_order = 3;
kernel_size = 11;
y_sg = sgolayfilt(y,polynomial_order,kernel_size);
y_sg_abs = abs(y_sg- offset);
[gfit] = fitGaussianNoOffset1D(x,y_sg_abs);
amp = (max(y_sg)-min(y_sg))/2;
[~,center_index] = min(abs(x));
if y_sg(center_index) > offset
    amp = -amp;
end
sig = gfit.sig;

% compute initial guesses

% define fit model
cosine = @(x,amp,k,phase) -amp*cos(k.*x+phase);
gaussian = @(x,cen,sig) exp(-(x-cen).^2./(2*sig^2));
fit_model = @(p,x) cosine(x,p(1),p(2),p(3)).*gaussian(x,p(4),p(5)) + p(6);

% initialize parameter constraints
p0 = []; lb = []; ub = [];
p0(1) = amp; lb(1) = -abs(amp)*10; ub(1) = abs(amp)*10;
p0(2) = k; lb(2) = k/2; ub(2) = k*2;
p0(3) = 0; lb(3) = -pi/4; ub(3) = pi/4;
p0(4) = 0; lb(4) = -sig/4; ub(4) = sig/4;
p0(5) = sig; lb(5) = sig/2; ub(5) = sig*2;
if nargin < 4
    p0(6) = 0; lb(6) = 0; ub(6) = 0;
else
    p0(6) = mean(y); lb(6) = min(y); ub(6) = max(y);
end


% do fitting and error estimation
p = lsqcurvefit(fit_model,p0,x,y,lb,ub);
% pci = nlparci(p,R,'jacobian',J);    % 95% confidence intervals from fit
% pse = (pci(:,2) - pci(:,1))/3.92;    % convert 95% confidence intervals to standard error

% output results
results.fit = fit_model(p,x);
results.guess = fit_model(p0,x);
results.amp = p(1);
results.k = p(2);
results.phase = p(3);
results.g_cen = p(4);
results.g_sig = p(5);
results.lam = 2*pi/results.k;
results.offset = p(6);

end