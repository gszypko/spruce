function [fit_p,Rx,Ry] = extractHoleOrientation(x,y,n,vis,theta,savedir)
% x (vector): position along dimension 2 of 'n' in cm
% y (vector): position along dimension 1 of 'n' in cm
% n (matrix): ion density in cm^-3
% vis (string): 'on' or 'off', controls visibility of figure
% savedir (string): full path to where figure is saved, not saved if empty

% default arguments
if nargin < 4, theta = 0; end
if nargin < 5, savedir = []; end

% 2D rotation transformation
Rx = @(x,y,t) x.*cosd(t) - y.*sind(t);
Ry = @(x,y,t) x.*sind(t) + y.*cosd(t);

% fit with ordinary 2D Gaussian to get initial guesses for overall density distribution
fit_g = fitImgWithGaussian(x,y,n,0);

% define fit model (2D Gaussian and hole with arbitrary principle axes
G2D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(Rx(x,y,t)-x0).^2./(2*sigx^2)-(Ry(x,y,t)-y0).^2./(2*sigy^2));
G1D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(Rx(x,y,t)-x0).^2./(2*sigx^2)-(Ry(x,y,t)-y0).^2./(2*sigy^2));
fun = @(c,d) G2D(d(:,1),d(:,2),c(1),c(4),c(2),c(3),c(5),c(6)) - G1D(d(:,1),d(:,2),c(9),c(7),c(2),c(3),c(8),c(6));

% initial guesses for fit model - derived from 2D Gaussian fit to density distribution
c0 = [theta fit_g.x0 fit_g.y0 fit_g.amp fit_g.sigx fit_g.sigy fit_g.amp/10 fit_g.sigx/5];
c0(1) = theta; lb(1) = c0(1)-30; ub(1) = c0(1)+30;
c0(2) = Rx(fit_g.x0,fit_g.y0,theta); lb(2) = c0(2)-.1; ub(2) = c0(2)+.1;
c0(3) = Ry(fit_g.x0,fit_g.y0,theta); lb(3) = c0(3)-.1; ub(3) = c0(3)+.1;
c0(4) = fit_g.amp; lb(4) = c0(4)/2; ub(4) = c0(4)*2;
c0(5) = fit_g.sigx; lb(5) = c0(5)/2; ub(5) = c0(5)*2;
c0(6) = fit_g.sigy; lb(6) = c0(6)/2; ub(6) = c0(6)*2;
c0(7) = fit_g.amp/10; lb(7) = c0(7)/5; ub(7) = c0(7)*5;
c0(8) = fit_g.sigx/8; lb(8) = c0(8)/5; ub(8) = c0(8)*5;
c0(9) = theta; lb(9) = c0(9)-30; ub(9) = c0(9)+30;


% prepare data for fit
ind_x = abs(x) < (2*fit_g.sigx - fit_g.x0);
ind_y = abs(y) < (2*fit_g.sigy - fit_g.y0);
[X,Y] = meshgrid(x(ind_x),y(ind_y));
xdata = [X(:) Y(:)];
ydata = n(ind_y,ind_x);
ydata = ydata(:);

% do fit
c = lsqcurvefit(fun,c0,xdata,ydata,lb,ub);

% unfold fit parameters
fit_p.amp_g = c(4);
fit_p.sigx_g = c(5);
fit_p.sigy_g = c(6);
fit_p.theta = c(9);
fit_p.amp = c(7);
fit_p.sigx = c(8);
fit_p.x0 = c(2);
fit_p.y0 = c(3);

% format fit data for plotting
p.x = x;
p.y = y;
[X,Y] = meshgrid(p.x,p.y);
xdata = [X(:) Y(:)];
p.n = reshape(n,length(p.y),length(p.x));
p.n0 = reshape(fun(c0,xdata),length(p.y),length(p.x));
p.nFit = reshape(fun(c,xdata),length(p.y),length(p.x));
p.nRes = sgfilt2D(p.nFit-p.n,5,5,3,3);

% plot fit data
q = {'n','n0','nFit','nRes'};
num = length(q);
[fig,ax,an,row,col] = open_subplot(num,'Visible',vis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
cmax = max(p.nFit,[],'all');
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        cax = get_axis(fig,ax{i,j});

        xdata = p.x;
        ydata = p.y;
        zdata = p.(q{iter});
        imagesc(xdata,ydata,zdata);
        colorbar
        hold on
        y = linspace(-3,3,101);
        plot(Rx(0,y,-fit_p.theta),Ry(0,y,-fit_p.theta))
        cax.YDir = 'normal';
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 12;
        if strcmp(q{iter},'nRes'), cax.CLim = [-1 1]*cmax/10;
        else, cax.CLim = [0 1]*cmax; end
        xlabel('x (cm)')
        ylabel('y (cm)')
        title(q{iter},'FontWeight','normal')
    end
end
an.String = ['Fit to Initial Time Point: \theta = ' num2str(fit_p.theta)];
an.Position = [0.1414    0.9399    0.7230    0.0592];
f = filesep;
if ~isempty(savedir), saveas(fig,[savedir f 'ion-hole-init.png']); end
close(fig)

end