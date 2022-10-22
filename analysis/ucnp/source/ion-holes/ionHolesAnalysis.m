function [] = ionHolesAnalysis(data,flags)
f = filesep;
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

%% Fit Gaussian to Initial Time Point to Get Hole Location

q = {'n','n0','nFit','nRes'};
num = length(q);
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
frames = cell(length(data.grids.time));
s(length(data.grids.time)).t = [];
for k = 1:length(data.grids.time)
    disp(['Fitting Ion Holes: ' num2str(k) '/' num2str(length(data.grids.time))])
    fit = fitImgWithGaussian(data.grids.x_vec,data.grids.y_vec,data.grids.vars(k).n);
    s(k).t = data.grids.time(k);

    x_prime = @(x,y,t) x.*cosd(t) - y.*sind(t);
    y_prime = @(x,y,t) x.*sind(t) + y.*cosd(t);
    G2D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(x_prime(x,y,t)-x0).^2./(2*sigx^2)-(y_prime(x,y,t)-y0).^2./(2*sigy^2));
    G1D = @(x,y,t,A,x0,y0,sigx,sigy) A.*exp(-(x_prime(x,y,t)-x0).^2./(2*sigx^2)-(y_prime(x,y,t)-y0).^2./(2*sigy^2));
    fun = @(c,d) G2D(d(:,1),d(:,2),c(1),c(4),c(2),c(3),c(5),c(6)) - G1D(d(:,1),d(:,2),c(1),c(7),c(2),c(3),c(8),c(6)) ...
                                                                  - G1D(d(:,1),d(:,2),c(1),c(7),c(9),c(3),c(8),c(6)) ...
                                                                  - G1D(d(:,1),d(:,2),c(1),c(7),c(10),c(3),c(8),c(6));
    
    % angle, xcen, ycen, amp1, sigx1, sigy1, amp2, sigx1
    if k == 1, s(1).theta = 15; end
    xp0 = x_prime(fit.x0,fit.y0,s(1).theta);
    sig0 = data.sig0;
    tau = data.tau;
    t_fac = sqrt(1+s(k).t^2/tau^2);
    xL = +(xp0+sig0*atan(s(k).t/tau))*t_fac;
    xR = -(xp0+sig0*atan(s(k).t/tau))*t_fac;
    x0 = [s(1).theta fit.x0 fit.y0 fit.amp fit.sigx fit.sigy fit.amp/10/3 fit.sigx/8 xL xR];
    lb(1) = s(1).theta-2; ub(1) = s(1).theta+2;
    lb(2) = -.1; ub(2) = .1;
    lb(3) = -.1; ub(3) = .1;
    lb(4) = fit.amp/2; ub(4) = fit.amp*2;
    lb(5) = fit.sigx/2; ub(5) = fit.sigx*2;
    lb(6) = fit.sigy/2; ub(6) = fit.sigy*2;
    lb(7) = fit.amp/100; ub(7) = fit.amp/5;
    lb(8) = fit.sigx/25; ub(8) = fit.sigx/4;
    lb(9) = xL-fit.sigx/2; ub(9) = xL+fit.sigx/2;
    lb(10) = xR-fit.sigx/2; ub(10) = xR+fit.sigx/2;
    if k~= 1
        x0(1) = s(1).theta;
        lb(1) = s(1).theta;
        ub(1) = s(1).theta;
    end
    
    [X,Y] = meshgrid(data.grids.x_vec,data.grids.y_vec);
    xdata = [X(:) Y(:)];
    [ydata] = sgfilt2D(data.grids.vars(k).n(:),15,15,3,3);
    
    x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
    
    s(k).theta = x(1);
    s(k).x0 = x(2);
    s(k).y0 = x(3);
    s(k).nmax = x(4);
    s(k).sigx = x(5);
    s(k).sigy = x(6);
    s(k).pmax = x(7);
    s(k).psig = x(8);
    s(k).xL = x(9);
    s(k).xR = x(10);
    
    s(k).x = data.grids.x_vec;
    s(k).y = data.grids.y_vec;
    s(k).n = reshape(ydata,length(s(k).x),length(s(k).y));
    s(k).n0 = reshape(fun(x0,xdata),length(s(k).x),length(s(k).y));
    s(k).nFit = reshape(fun(x,xdata),length(s(k).x),length(s(k).y));
    s(k).nRes = s(k).nFit - s(k).n;
    
    
    cmax = [min(s(k).nFit,[],'all') max(s(k).nFit,[],'all')];
    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            cax = get_axis(fig,ax{i,j});
    
            xdata = s(k).x;
            ydata = s(k).y;
            zdata = s(k).(q{iter});
            imagesc(xdata,ydata,zdata);
            colorbar
    
            cax.YDir = 'normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 12;
            cax.CLim = cmax;
            if strcmp(q{iter},'nRes'), cax.CLim = cmax/10; end
            xlabel('x (cm)')
            ylabel('y (cm)')
            title(q{iter},'FontWeight','normal')
        end
    end
    dlm = ' - ';
    str1 = ['Iter = ' num2str(k-1)];
    str2 = ['t = ' num2str(data.grids.vars(k).time*1e6,'%.3g') '\mus = ' num2str(data.grids.vars(k).time/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = [str1 dlm str2];
    an.Position = [0.1636    0.9308    0.7230    0.0801];
    frames{k} = getframe(fig);
end
close(fig)
write_video([data.folder f 'ion-hole-fits-2D-n'],frames);

for k = 1:length([s.t])
    t_fac = sqrt(1+s(k).t^2/data.tau^2);
    xp0 = x_prime(fit.x0,fit.y0,s(1).theta);
    s(k).xL_model = -(xp0+data.sig0*atan(s(k).t/tau))*t_fac;
    s(k).xR_model = +(xp0+data.sig0*atan(s(k).t/tau))*t_fac;
end

num = 1;
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
q = {'xL','xR','xL_model','xR_model'};
qstr = {'xL','xR','Model L','Model R'};
l = get_line_specs(length(q));
cax = get_axis(fig,ax{1});
hold on
for i = 1:length(q)
    plot([s.t]/data.tau,[s.(q{i})],'LineWidth',2,'MarkerSize',4,'Color',l(i).col,'MarkerFaceColor',l(i).col,'MarkerEdgeColor',l(i).col)
end
cax.FontSize = 12;
xlabel('t / \tau')
ylabel('Hole Pos. (cm)')
title('Extracted from 2D Fits to Ion Density')
lgd = legend(qstr);
lgd.Position = [0.1865    0.7103    0.2306    0.1921];
saveas(fig,[data.folder f 'ion-hole-pos-2D-n.png']);
close(fig)


%% 


end