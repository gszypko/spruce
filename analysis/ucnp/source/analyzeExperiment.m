function [gauss_fit_int,vel_data,gauss_fit_rot] = analyzeExperiment(path,flags)
%% Load Experimental Data
f = filesep;
load([path f 'os.mat'],'os');
[~,ind_t] = sort([os.delays]); 
os = os(ind_t);
for i = 1:length(os)
    os(i).delays = os(i).delays+os(i).tE/2;
end

%% Obtain Density from Integrated Spectra
imgs_int(length([os.delays])).t = [];
for i = 1:length([os.delays])
    % retrieve information from output structure
    imgs_int(i).t = [os(i).delays]*1e-9;
    imgs_int(i).x = os(i).imgs.xRelInMM./10;
    imgs_int(i).y = os(i).imgs.yRelInMM./10;
    imgs_int(i).n = os(i).imgs.density.*1e8;

    % trim data limits
    ind_x = abs(imgs_int(i).x) < flags.trim_exp_domain(1);
    ind_y = abs(imgs_int(i).y) < flags.trim_exp_domain(2);
    imgs_int(i).x = imgs_int(i).x(ind_x);
    imgs_int(i).y = imgs_int(i).y(ind_y);
    imgs_int(i).n = imgs_int(i).n(ind_y,ind_x);
end

%% Extract Plasma Size from Density Distribution

% extract sizes from fit
gauss_fit_int = struct;
gauss_fit_int(length([imgs_int.t])).sig_x = [];
gauss_fit_int(length([imgs_int.t])).sig_y = [];
for i = 1:length([imgs_int.t])
    results = fitImgWithGaussian(imgs_int(i).x,imgs_int(i).y,imgs_int(i).n,0);
    gauss_fit_int(i).t = imgs_int(i).t;
    gauss_fit_int(i).x = results.x;
    gauss_fit_int(i).y = results.y;
    gauss_fit_int(i).img = results.img;
    gauss_fit_int(i).imgfit = results.imgfit;
    gauss_fit_int(i).imgres = results.imgres;
    gauss_fit_int(i).sig_x = results.sigx;
    gauss_fit_int(i).sig_y = results.sigy;
    gauss_fit_int(i).sig_2D = (gauss_fit_int(i).sig_x*gauss_fit_int(i).sig_y)^(1/2); % 1/3 because this is the experimental data
    gauss_fit_int(i).sig_3D = (gauss_fit_int(i).sig_x*gauss_fit_int(i).sig_y^2)^(1/3); % 1/3 because this is the experimental data
end

%% Plot Gaussian Fits to Experimental Density
vars = {'imgfit','img','imgres'};
num = length(vars);
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.159500000000000,0.887715907879935,0.723000000000000,0.080100000000000];
for k = 1:length(gauss_fit_int) 
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end

            cax = get_axis(fig,ax{i,j});
            xdata = gauss_fit_int(k).x;
            ydata = gauss_fit_int(k).y;
            zdata = gauss_fit_int(k).(vars{iter});
            imagesc(xdata,ydata,zdata)
            colorbar
            cax.YDir = 'Normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            if i == size(ax,1), xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(vars{iter},'FontWeight','normal')
            if iter == 1, clim = max(zdata,[],'all'); end
            ax{i,j}.CLim = [0 clim];
            if iter == 3, ax{i,j}.CLim = [-clim clim]./5; end
        end
    end
    
    str1 = 'Gaussian Fits to Images: ';
    str2 = ['t = ' num2str(gauss_fit_int(k).t*1e6,'%.3g') '\mus'];
    an.String = [str1 str2];
    saveas(fig,[path f 'gauss-fit-' num2str(k) '.png']);
end
close(fig)


%% Extract Te from Size Evolution

% extract Te(0) from size evolution along x axis
tau = @(sig,T) getTauExp(sig,T);
fitmodel = @(c,data) c(1)*sqrt(1+data.^2/tau(c(1),c(2))^2);
p0 = [gauss_fit_int.sig_x os(1).Te];
data = [imgs_int.t];
data = data';
zdata = [gauss_fit_int.sig_x]';
[cx] = lsqcurvefit(fitmodel,p0,data,zdata,0,200);
Te_from_sig_x = cx(2);

% extract Te(0) from size evolution along x axis
p0 = [gauss_fit_int.sig_y os(1).Te];
data = [imgs_int.t];
data = data';
zdata = [gauss_fit_int.sig_y]';
[cy] = lsqcurvefit(fitmodel,p0,data,zdata,0,200);
Te_from_sig_y = cy(2);

% extract Te(0) from geometric mean of plasma size
p0 = [gauss_fit_int(1).sig_2D os(1).Te];
data = [imgs_int.t];
data = data';
zdata = [gauss_fit_int.sig_2D]';
[c2D] = lsqcurvefit(fitmodel,p0,data,zdata,0,200);
Te_from_sig_2D = c2D(2);

% open figure
num = 3;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
fig.Position = [22.200000000000003,2.666000000000000e+02,9.376000000000000e+02,2.951999999999999e+02];

% plot x axis information
[cax] = get_axis(fig,ax{1}); %#ok<NASGU> 
l = get_line_specs(2);  
hold on
plot(data.*1e6,fitmodel(cx,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data.*1e6,[gauss_fit_int.sig_x],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
xlabel('t (\mus)')
ylabel('\sigma_x (cm)')
title(['T_e = ' num2str(Te_from_sig_x) ' K'],'FontWeight','normal')
cax.PlotBoxAspectRatio = [1 1 1];
grid on 
grid minor

% plot y axis information
[cax] = get_axis(fig,ax{2}); %#ok<NASGU> 
l = get_line_specs(2);  
hold on
plot(data.*1e6,fitmodel(cy,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data.*1e6,[gauss_fit_int.sig_y],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
xlabel('t (\mus)')
ylabel('\sigma_y (cm)')
cax.PlotBoxAspectRatio = [1 1 1];
title(['T_e = ' num2str(Te_from_sig_y) ' K'],'FontWeight','normal')
grid on 
grid minor
lgd = legend({'data','fit'});
lgd.Position = [0.814002811768537,0.137802854962636,0.078427250381274,0.098215042129087];

% plot y axis information
[cax] = get_axis(fig,ax{3}); %#ok<NASGU> 
l = get_line_specs(2);  
hold on
plot(data.*1e6,fitmodel(c2D,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data.*1e6,[gauss_fit_int.sig_2D],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
xlabel('t (\mus)')
ylabel('\sigma (cm)')
title(['T_e = ' num2str(Te_from_sig_2D) ' K'],'FontWeight','normal')
cax.PlotBoxAspectRatio = [1 1 1];
grid on 
grid minor
lgd = legend({'data','fit'});
lgd.Position = [0.820828750335090,0.219938830148715,0.078427250381274,0.112804880374815];
an.String = 'Evolution of RMS Plasma Size';
saveas(fig,[path f 'Te_cartesian.png'])
close(fig)

%% Extract Te(0) from Velocity Transects

% get data from output structure
vel_data = struct();
vel_data(length([imgs_int.t])).t = [];
vel_data(length([imgs_int.t])).x = [];
vel_data(length([imgs_int.t])).v = [];
for i = 1:length(vel_data)
    vel_data(i).t = os(i).delays*1e-9;
    vel_data(i).x = [os(i).local.x]./10;
    vel_data(i).v = [os(i).local.vExp]*100;
end

% fit velocity distribution along x axis for each time point to extract Te(0)
tau_x = @(T) tau(gauss_fit_int(1).sig_x,T);
gam_x = @(t,T) (t/tau_x(T)^2)/(1+t^2/tau_x(T)^2);
for i = 1:length(vel_data)
    fitmodel = @(c,data) data.*gam_x(vel_data(i).t,c);
    c0 = Te_from_sig_x;
    ind = abs(vel_data(i).x) < 0.3;
    data = [vel_data(i).x(ind)];
    zdata = [vel_data(i).v(ind)];
    vel_data(i).Te_fit = lsqcurvefit(fitmodel,c0,data,zdata,1,200);
    vel_data(i).v_fit = fitmodel(vel_data(i).Te_fit,vel_data(i).x);
end

% opening a subplot 
num = length(vel_data);
[fig,ax,an] = open_subplot(num-1);

% updating axis information
vars = {'v','v_fit'};
iter = 1;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end
        
        cax = get_axis(fig,ax{i,j});
        l = get_line_specs(length(vars));
        for k = 1:length(vars)
            xdata = [vel_data(iter).x];
            ydata = [vel_data(iter).(vars{k})];
            plot(xdata,ydata,'o','LineWidth',2,'MarkerSize',2,'Color',l(k).col,'MarkerFaceColor',l(k).col,'MarkerEdgeColor',l(k).col)
            hold on
        end
        cax.FontSize = 11;
        if i == size(ax,1), xlabel('x (cm)'), end
        if j == 1, ylabel('v (cm/s)'), end
        title(['t = ' num2str(vel_data(iter).t*1e6) ' \mus'],'FontWeight','normal')
    end
end
Te = [vel_data.Te_fit];
Te_avg = mean(Te(2:end));
an.String = ['Mean Te from Fit: ' num2str(Te_avg)];
lgd = legend(vars);
lgd.Position = [0.8117    0.1192    0.0658    0.0894];
saveas(fig,[path f 'velocity-transects.png'])
close(fig)

%% Determine Ion Hole Orientation
[~,ind_t] = min([imgs_int.t]);
x = imgs_int(ind_t).x;
y = imgs_int(ind_t).y;
n = imgs_int(ind_t).n;
[hole_fit,Rx,Ry] = extractHoleOrientation(x,y,n,flags.figvis,0,path);

%% Project Data on Grid Aligned with Ion Hole Propagation Axis
% Uses grid transformation from 'extractHoleOrientation'
imgs_rot = struct;
imgs_rot(length([imgs_int.t])).t = [];
for k = 1:length([imgs_int.t])
    % record current time
    imgs_rot(k).t = imgs_int(k).t;

    % original spatial grids
    x = imgs_int.x;
    y = imgs_int.y;
    [X,Y] = meshgrid(x,y);
    
    % rotated spatial grids
    x_rot = imgs_int(k).x-hole_fit.x0;
    y_rot = imgs_int(k).y-hole_fit.y0;
    [X_rot,Y_rot] = meshgrid(x_rot,y_rot);
    
    % rotate grids to obtain values along the ion hole propagation axis - perpendicular quantities are ignored
    imgs_rot(k).x = x_rot;
    imgs_rot(k).y = y_rot;
    imgs_rot(k).n = interp2(X,Y,imgs_int(k).n,Rx(X_rot,Y_rot,-hole_fit.theta),Ry(X_rot,Y_rot,-hole_fit.theta),'linear',0);    
end

%% Extract Plasma Size from Rotated Density Distribution

% extract sizes from fit
gauss_fit_rot = struct;
gauss_fit_rot(length([imgs_rot.t])).sig_x = [];
gauss_fit_rot(length([imgs_rot.t])).sig_y = [];
for i = 1:length([imgs_rot.t])
    results = fitImgWithGaussian(imgs_rot(i).x,imgs_rot(i).y,imgs_rot(i).n,0);
    gauss_fit_rot(i).sig_x = results.sigx;
    gauss_fit_rot(i).sig_y = results.sigy;
    gauss_fit_rot(i).sig_2D = (gauss_fit_rot(i).sig_x*gauss_fit_rot(i).sig_y)^(1/2); % 1/3 because this is the experimental data
    gauss_fit_rot(i).sig_3D = (gauss_fit_rot(i).sig_x*gauss_fit_rot(i).sig_y^2)^(1/3); % 1/3 because this is the experimental data
end

%% Extract Te from Size Evolution

% extract Te(0) from size evolution along x axis
tau_x = @(T) getTauExp(gauss_fit_rot(1).sig_x,T);
fitmodelx = @(c,data) gauss_fit_rot(1).sig_x*sqrt(1+data.^2/tau_x(c)^2);
p0 = os(1).Te;
data = [imgs_int.t];
data = data';
zdata = [gauss_fit_rot.sig_x]';
Te_from_sig_x = lsqcurvefit(fitmodelx,p0,data,zdata,0,200);

% extract Te(0) from size evolution along x axis
tau_y = @(T) getTauExp(gauss_fit_rot(1).sig_y,T);
fitmodely = @(c,data) gauss_fit_rot(1).sig_y*sqrt(1+data.^2/tau_y(c)^2);
p0 = os(1).Te;
data = [imgs_int.t];
data = data';
zdata = [gauss_fit_rot.sig_y]';
Te_from_sig_y = lsqcurvefit(fitmodely,p0,data,zdata,0,200);

% open figure
num = 2;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
fig.Position = [2.142000000000000e+02,2.258000000000000e+02,8.724000000000001e+02,4.208000000000000e+02];
an.Position = [0.159500000000000,0.951400000000000,0.723000000000000,0.037100000000000];

% plot x axis information
[cax] = get_axis(fig,ax{1}); %#ok<NASGU> 
l = get_line_specs(2);  
hold on
plot(data.*1e6,fitmodelx(Te_from_sig_x,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data.*1e6,[gauss_fit_rot.sig_x],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
xlabel('t (\mus)')
ylabel('\sigma_x (cm)')
title(['T_e = ' num2str(Te_from_sig_x) ' K'],'FontWeight','normal')
grid on 
grid minor

% plot y axis information
[cax] = get_axis(fig,ax{2}); %#ok<NASGU> 
l = get_line_specs(2);  
hold on
plot(data.*1e6,fitmodely(Te_from_sig_y,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data.*1e6,[gauss_fit_rot.sig_y],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
xlabel('t (\mus)')
ylabel('\sigma_y (cm)')
title(['T_e = ' num2str(Te_from_sig_y) ' K'],'FontWeight','normal')
grid on 
grid minor
lgd = legend({'data','fit'});
lgd.Position = [0.807709788310235,0.192706836263204,0.085511234024149,0.083464113748483];
an.String = 'Evolution of RMS Plasma Radius Along Principal Axes';
saveas(fig,[path f 'Te_principal.png'])
close(fig)

end