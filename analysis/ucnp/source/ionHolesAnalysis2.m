function [] = ionHolesAnalysis2(data,flags)
f = filesep;
theta = 20; % initial guess for rotation angle that defines principle axes (degree)

%% Check for Required Grids
% In this section, I compute any grids that I would like to work with that were not saved in the .out file
% I am expecting for grids n, v_x, v_y, i_temp, e_temp, dt, dPdx, dPdy to exist already, so I will run a check for that
expected_grids = {'n','v_x','v_y','i_temp','e_temp','dt','dPdx','dPdy'};
vars = fieldnames(data.grids.vars);
for i = 1:length(expected_grids)
    if ~max(strcmp(vars,expected_grids{i}))
        error(['Expected variable <' expected_grids{i} '> is missing.'])
    end
end

%% Fit Gaussian with Hole to Initial Time Point

% fix rotated coordinate system
x_rot = @(x,y,t) x.*cosd(t) - y.*sind(t);
y_rot = @(x,y,t) x.*sind(t) + y.*cosd(t);

% obtain index to smallest time point
[~,ind_t] = min(data.grids.time);
density = data.grids.vars(1).n;

% fit with ordinary 2D Gaussian to get initial guesses for overall density distribution
fit = fitImgWithGaussian(data.grids.x_vec,data.grids.y_vec,density,true);

% define fit model (2D Gaussian and hole with arbitrary principle axes
G2D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(x_rot(x,y,t)-x0).^2./(2*sigx^2)-(y_rot(x,y,t)-y0).^2./(2*sigy^2));
G1D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(x_rot(x,y,t)-x0).^2./(2*sigx^2)-(y_rot(x,y,t)-y0).^2./(2*sigy^2));
fun = @(c,d) G2D(d(:,1),d(:,2),c(1),c(4),c(2),c(3),c(5),c(6)) - G1D(d(:,1),d(:,2),c(1),c(7),c(2),c(3),c(8),c(6));

% initial guesses for fit model - derived from 2D Gaussian fit to density distribution
c0 = [theta fit.x0 fit.y0 fit.amp fit.sigx fit.sigy fit.amp/10 fit.sigx/5];
c0(1) = theta; lb(1) = theta-10; ub(1) = theta+10;
c0(2) = x_rot(fit.x0,fit.y0,theta); lb(2) = c0(2)-.1; ub(2) = c0(2)+.1;
c0(3) = y_rot(fit.x0,fit.y0,theta); lb(3) = c0(3)-.1; ub(3) = c0(3)+.1;
c0(4) = fit.amp; lb(4) = fit.amp/1.25; ub(4) = fit.amp*1.25;
c0(5) = fit.sigx; lb(5) = fit.sigx/1.25; ub(5) = fit.sigx*1.25;
c0(6) = fit.sigy; lb(6) = fit.sigy/1.25; ub(6) = fit.sigy*1.25;
c0(7) = fit.amp/10; lb(7) = fit.amp/10/3; ub(7) = fit.amp/10*3;
c0(8) = fit.sigx/8; lb(8) = fit.sigx/8/3; ub(8) = fit.sigx/8*3;

% prepare data for fit
ind_x = abs(data.grids.x_vec) < 2*fit.sigx;
ind_y = abs(data.grids.y_vec) < 2*fit.sigy;
[X,Y] = meshgrid(data.grids.x_vec(ind_x),data.grids.y_vec(ind_y));
xdata = [X(:) Y(:)];
ydata = density(ind_y,ind_x);
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
p.n = reshape(density,length(p.x),length(p.y));
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

%% Project Data on Grid Aligned with Perturbation
prop = struct;
prop(length(data.grids.time)).t = [];
num = 1;
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
fig.Position = [257.0000  234.6000  784.0000  641.6000];
for k = 1:length(data.grids.time)
    % record current time
    prop(k).t = data.grids.time(k);

    % project density grid onto rotated coordinate system
    xp = data.grids.x_vec-hole.xp0;
    yp = data.grids.y_vec-hole.yp0;
    [Xp,Yp] = meshgrid(xp,yp);
    np = interp2(X,Y,data.grids.vars(k).n,x_rot(Xp,Yp,-hole.theta),y_rot(Xp,Yp,-hole.theta),'spline',min(data.grids.vars(k).n,[],'all'));

    % fit with 2D Gaussian
    [fit] = fitImgWithGaussian(xp,yp,np,false);
    ind_y = abs(yp) < fit.sigy/2;
    yp = yp(ind_y);
    np = np(ind_y,:);
    np = mean(np,1);
    

    % plot density distribution along rotated axis
    
    plot(xp,np)
    pause(0.5)
end



end