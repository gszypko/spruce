function [fit,x_rot,y_rot,x_rot_inv,y_rot_inv] = extractHoleOrientation(x,y,n)
% x (vector): position along dimension 2 of 'n' in cm
% y (vector): position along dimension 1 of 'n' in cm
% n (matrix): plasma/ion density in cm^-^3

% 2D rotation transformation
x_rot = @(x,y,t) x.*cosd(t) - y.*sind(t);
y_rot = @(x,y,t) x.*sind(t) + y.*cosd(t);

% fit with ordinary 2D Gaussian to get initial guesses for overall density distribution
fit = fitImgWithGaussian(x,y,n,true);

% define fit model (2D Gaussian and hole with arbitrary principle axes
G2D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(x_rot(x,y,t)-x0).^2./(2*sigx^2)-(y_rot(x,y,t)-y0).^2./(2*sigy^2));
G1D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(x_rot(x,y,t)-x0).^2./(2*sigx^2)-(y_rot(x,y,t)-y0).^2./(2*sigy^2));
fun = @(c,d) G2D(d(:,1),d(:,2),c(1),c(4),c(2),c(3),c(5),c(6)) - G1D(d(:,1),d(:,2),c(1),c(7),c(2),c(3),c(8),c(6));

% initial guesses for fit model - derived from 2D Gaussian fit to density distribution
c0 = [theta fit.x0 fit.y0 fit.amp fit.sigx fit.sigy fit.amp/10 fit.sigx/5];
c0(1) = theta; lb(1) = c0(1)-10; ub(1) = c0(1)+10;
c0(2) = x_rot(fit.x0,fit.y0,theta); lb(2) = c0(2)-.1; ub(2) = c0(2)+.1;
c0(3) = y_rot(fit.x0,fit.y0,theta); lb(3) = c0(3)-.1; ub(3) = c0(3)+.1;
c0(4) = fit.amp; lb(4) = c0(4)/1.25; ub(4) = c0(4)*1.25;
c0(5) = fit.sigx; lb(5) = c0(5)/1.25; ub(5) = c0(5)*1.25;
c0(6) = fit.sigy; lb(6) = c0(6)/1.25; ub(6) = c0(6)*1.25;
c0(7) = c0(7)/10; lb(7) = c0(7)/10/3; ub(7) = c0(7)/10*3;
c0(8) = fit.sigx/8; lb(8) = c0(8)/8/3; ub(8) = c0(8)/8*3;

% prepare data for fit
ind_x = abs(data.grids.x_vec) < 2*fit.sigx;
ind_y = abs(data.grids.y_vec) < 2*fit.sigy;
[X,Y] = meshgrid(data.grids.x_vec(ind_x),data.grids.y_vec(ind_y));
xdata = [X(:) Y(:)];
ydata = n(ind_y,ind_x);
ydata = ydata(:);

% do fit
c = lsqcurvefit(fun,c0,xdata,ydata,lb,ub);
hole.theta = c(1);
hole.xp0 = c(2);
hole.yp0 = c(3);
hole.amp = c(7);
hole.sigx = c(8);
hole.sigy = c(6);

% format fit data for plotting
p.x = data.grids.x_vec;
p.y = data.grids.y_vec;
[X,Y] = meshgrid(p.x,p.y);
xdata = [X(:) Y(:)];
p.n = reshape(n,length(p.x),length(p.y));
p.n0 = reshape(fun(c0,xdata),length(p.x),length(p.y));
p.nFit = reshape(fun(c,xdata),length(p.x),length(p.y));
p.nRes = p.nFit - p.n;

% plot fit data
q = {'n','n0','nFit','nRes'};
num = length(q);
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
fig.Position = [257.0000  234.6000  784.0000  641.6000];
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
        plot(x_rot(0,y,-hole.theta),y_rot(0,y,-hole.theta))
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
an.String = ['Fit to Initial Time Point: \theta = ' num2str(hole.theta)];
an.Position = [0.1414    0.9399    0.7230    0.0592];

end