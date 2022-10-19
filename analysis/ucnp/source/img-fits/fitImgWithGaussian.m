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
[xfilt,yfilt,imgfilt] = sgfilt2D(xbin,ybin,imgbin,x_sg,y_sg,3,3,true);
normfac = max(imgfilt,[],'all');
imgbin = imgbin./normfac;
imgfilt = imgfilt./normfac;

%% Fit Image with 2D Gaussian

% define fit model
% c: [n0 x0 y0 sigX sigY offset]
fitmodel = @(c,data) gaussian2D(c,data);

% initial guess for plasma center
imgint = sum(imgfilt,'all');
x0 = sum(sum(imgfilt,1).*xfilt,2)/imgint;
y0 = sum(sum(imgfilt,2).*yfilt',1)/imgint;

% initial guess for amplitude
offset = min(imgfilt,[],'all');
amp = max(imgfilt,[],'all') - offset;

% initial guess for widths
[X,Y] = meshgrid(xbin,ybin);
grid = [X(:) Y(:)];
ind = imgfilt(:) > (amp/2 + offset);
sigx = (max(grid(ind,1)) - min(grid(ind,1)))/(2*sqrt(2*log(2)));
sigy = (max(grid(ind,2)) - min(grid(ind,2)))/(2*sqrt(2*log(2)));

% define initial guesses and parameter bounds
p0(1) = amp; lb(1) = amp/5; ub(1) = amp*5;
p0(2) = x0; lb(2) = min(xbin); ub(2) = max(xbin);
p0(3) = y0; lb(3) = min(ybin); ub(3) = max(ybin);
p0(4) = sigx; lb(4) = sigx/5; ub(4) = sigx*5;
p0(5) = sigy; lb(5) = sigy/5; ub(5) = sigy*5;
p0(6) = offset; lb(6) = offset - amp/10; ub(6) = offset + amp/10; 
if offset_opt
    p0(6) = 0; lb(6) = 0; ub(6) = 0;
end

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
if offset_opt
    fit.offset = 0;
    fit.offsetErr = 0;
else 
    fit.offset = p(6)*normfac;
    fit.offsetErr = pse(6)*normfac;
end

end