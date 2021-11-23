function [fit] = fitImgWithGaussian(x,y,img)
% x (1xn double): x position in mm, varies with column number
% y (1xm double): y position in mm, varies with row number
% img (mxn double): image img(m,n) corresponds to x(n) and y(m)

%% Bin and Normalize Image For Speed
[xbin,ybin,imgbin] = binImgForFit(x,y,img,300);
imgfilt = imgaussfilt(imgbin,2);
normfac = max(imgfilt,[],'all');
imgbin = imgbin./normfac;
imgfilt = imgfilt./normfac;

%% Obtain Initial Guesses
% initial guess for plasma center
imgint = sum(imgfilt,'all');
x0 = sum(sum(imgfilt,1).*xbin,2)/imgint;
y0 = sum(sum(imgfilt,2).*ybin',1)/imgint;

% initial guess for amplitude
offset = min(imgfilt,[],'all');
amp = max(imgfilt,[],'all') - offset;

% initial guess for widths
[X,Y] = meshgrid(xbin,ybin);
grid = [X(:) Y(:)];
ind = imgfilt(:) > amp/2 + offset;
sigx = (max(grid(ind,1)) - min(grid(ind,1)))/(2*sqrt(2*log(2)));
sigy = (max(grid(ind,2)) - min(grid(ind,2)))/(2*sqrt(2*log(2)));

% define initial guesses and parameter bounds
p0 = [abs(amp) x0 y0 sigx sigy offset];

lb(1) = p0(1)/10; ub(1) = p0(1)*10;
lb(2) = min(xbin); ub(2) = max(xbin);
lb(3) = min(ybin); ub(3) = max(ybin);
lb(4) = p0(4)/5; ub(4) = p0(4)*5;
lb(5) = p0(5)/5; ub(5) = p0(5)*5;
% lb(6) = p0(6) - p0(1)/10; ub(6) = p0(6) + p0(1)/10; 


%% Do 2D Gaussian Fit
% Fit Model
% c: [n0 x0 y0 sigX sigY offset]
fitmodel = @(c,data) gaussian2D(c,data,0);

% Format data for fit
data = grid;
zdata = imgbin(:);

% Set fit options to default
options = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective');

% Do the fit
[p,~,R,~,~,~,J] = lsqcurvefit(fitmodel,p0,data,zdata,lb,ub,options);

pci = nlparci(p,R,'jacobian',J);    % 95% confidence intervals from fit
pse = (pci(:,2) - pci(:,1))/3.92;    % convert 95% confidence intervals to standard error

%% Output Fit Results
fit.xRelInMM = xbin;
fit.yRelInMM = ybin;
fit.imgbin = imgbin.*normfac;
fit.imgfit = reshape(fitmodel(p,data),size(imgbin,1),size(imgbin,2)).*normfac;
fit.imgres = fit.imgfit - fit.imgbin;

fit.amp = p(1).*normfac;
fit.ampErr = pse(1).*normfac;
fit.x0 = p(2);
fit.x0Err = pse(2);
fit.y0 = p(3);
fit.y0Err = pse(3);
fit.sigx = p(4);
fit.sigxErr = pse(4);
fit.sigy = p(5);
fit.sigyErr = pse(5);
if nargin == 3
    fit.offset = offset.*normfac;
    fit.offsetErr = 0;
else 
    fit.offset = p(6);
    fit.offsetErr = pse(6);
end

end