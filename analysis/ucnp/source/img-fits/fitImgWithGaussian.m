function [fit] = fitImgWithGaussian(x,y,img,offset_opt)
% x (1xn double): x pos - varies with column number
% y (1xm double): y pos - varies with row number
% img (mxn double): image img(m,n) corresponds to x(n) and y(m)
% offset_opt (boolean): (true) fix offset to zero (false) fit offset

if nargin < 4, offset_opt = false; end

%% Bin and Normalize Image For Speed
[xbin,ybin,imgbin] = binImgForFit(x,y,img,250);
x_sg = ceil(size(img,2)/10);
y_sg = ceil(size(img,1)/10);
if x_sg < 5, x_sg = 5; end
if y_sg < 5, y_sg = 5; end
[imgfilt,ind_x,ind_y] = sgfilt2D(imgbin,x_sg,y_sg,3,3,true);
xfilt = xbin(ind_x);
yfilt = ybin(ind_y);
normfac = max(imgfilt,[],'all');
imgbin = imgbin./normfac;
imgfilt = imgfilt./normfac;

%% Fit Image with 2D Gaussian

% define fit model
% c: [n0 x0 y0 sigX sigY offset]
fitmodel = @(c,data) gaussian2D([c 0],data);

% initial guess for plasma center
imgintx = trapz(xfilt,imgfilt,2);
imginty = trapz(yfilt,imgfilt,1);
imgintxy = trapz(yfilt,trapz(xfilt,imgfilt,2));
x0 = trapz(xfilt,xfilt.*imginty)/imgintxy;
y0 = trapz(yfilt,yfilt.*imgintx')/imgintxy;

% initial guess for amplitude
offset = min(imgfilt,[],'all');
amp = max(imgfilt,[],'all') - offset;

% initial guess for widths
[X,Y] = meshgrid(xbin,ybin);
grid = [X(:) Y(:)];
ind = imgfilt(:) > (amp/2 + offset);
sigx = sqrt(trapz(xfilt,(xfilt-x0).^2.*imginty)/imgintxy);
sigy = sqrt(trapz(yfilt,(yfilt-y0).^2.*imgintx')/imgintxy);

% define initial guesses and parameter bounds
p0(1) = amp; lb(1) = amp/5; ub(1) = amp*5;
p0(2) = x0; lb(2) = min(xbin)/3; ub(2) = max(xbin)/3;
p0(3) = y0; lb(3) = min(ybin)/3; ub(3) = max(ybin)/3;
p0(4) = sigx; lb(4) = sigx/5; ub(4) = sigx*5;
p0(5) = sigy; lb(5) = sigy/5; ub(5) = sigy*5;

% format data for fit
data = grid;
zdata = imgbin(:);

% do the fit
[p,~,R,~,~,~,J] = lsqcurvefit(fitmodel,p0,data,zdata,lb,ub);

% obtain error statistics
pci = nlparci(p,R,'jacobian',J);    % 95% confidence intervals from fit
pse = (pci(:,2) - pci(:,1))/3.92;    % convert 95% confidence intervals to standard error

% output fit results
fit.x = x;
fit.y = y;
fit.img = img;
[X,Y] = meshgrid(fit.x,fit.y);
fit.imgguess = reshape(fitmodel(p0,[X(:) Y(:)]),size(img,1),size(img,2)).*normfac;
fit.imgfit = reshape(fitmodel(p,[X(:) Y(:)]),size(img,1),size(img,2)).*normfac;
fit.imgres = fit.imgfit - fit.img;

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

end