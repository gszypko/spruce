function [s] = compareExpAndSimData(data,flags)
%% load experimental data from <os.mat> file
% load data
f = filesep;
base_folder = extractBefore(data.folder,[f 'set_']);
load([base_folder f 'os.mat'],'os');
os = os([os.delays]*1e-9 < max(data.grids.time));

% find indices of simulation time points that are closest to experimental time points
% ind_t(i) holds the index of data.grids.vars corresponding to experimental time point os(i).delays
closest_time_point = zeros(size(os)); 
for i = 1:length([os.delays])
    [~,closest_time_point(i)] = min(abs(os(i).delays*1e-9 + os(i).tE/2*1e-9 - data.grids.time));
end

%% create structure for 2D image info and fit with 2D Gaussian
% obtain domain limits
x_lim = flags.plot_window(1);
y_lim = flags.plot_window(2);

% create struct to hold info to be plotted - all structures indexed by time
t = [os.delays].*1e-9;
t_fac = sqrt(1+t.^2/data.tau^2);
imgs = struct();
imgs(length([os.delays])).x = [];
tr = struct();
tr(length(imgs)).x = [];
sig = struct();
sig(length(imgs)).exp_x = [];

for i = 1:length([os.delays])    
    % trim grids for plotting
    x = os(i).imgs.xRelInMM/10;
    y = os(i).imgs.yRelInMM/10;
    ind_x = abs(x) < x_lim;
    ind_y = abs(y) < y_lim;
    imgs(i).x = x(ind_x);
    imgs(i).y = y(ind_y);
    [imgs(i).X, imgs(i).Y] = meshgrid(imgs(i).x, imgs(i).y);
    imgs(i).n_exp = os(i).imgs.density(ind_y,ind_x)*1e8;
    imgs(i).n_sim = interp2(data.grids.pos_x,data.grids.pos_y,data.grids.vars(closest_time_point(i)).n,imgs(i).X,imgs(i).Y,'linear',0);
    imgs(i).n_res = imgs(i).n_sim - imgs(i).n_exp;
    imgs(i).v_sim_x = interp2(data.grids.pos_x,data.grids.pos_y,data.grids.vars(closest_time_point(i)).v_x,imgs(i).X,imgs(i).Y,'linear',0);
    imgs(i).v_sim_y = interp2(data.grids.pos_x,data.grids.pos_y,data.grids.vars(closest_time_point(i)).v_y,imgs(i).X,imgs(i).Y,'linear',0);
    
    % get plasma sizes
    imgs(i).fit_exp = fitImgWithGaussian(imgs(i).x,imgs(i).y,imgs(i).n_exp,0);
    imgs(i).fit_sim = fitImgWithGaussian(imgs(i).x,imgs(i).y,imgs(i).n_sim,0);
    sig(i).exp_x = imgs(i).fit_exp.sigx;
    sig(i).sim_x = imgs(i).fit_sim.sigx;
    sig(i).exp_y = imgs(i).fit_exp.sigy;
    sig(i).sim_y = imgs(i).fit_sim.sigy;
    sig(i).exp_2D = sqrt(imgs(i).fit_exp.sigx*imgs(i).fit_exp.sigy);
    sig(i).sim_2D = sqrt(imgs(i).fit_sim.sigx*imgs(i).fit_sim.sigy);
    sig(i).exp_3D = (imgs(i).fit_exp.sigx*imgs(i).fit_exp.sigy^2)^(1/3);

    % get transect data
    [~,ind_y] = min(abs(imgs(i).y));
    [~,ind_x] = min(abs(imgs(i).x));
    tr(i).x = imgs(i).x;
    tr(i).y = imgs(i).y;

    tr(i).n_exp = imgs(i).n_exp(ind_y,:);
    tr(i).n_sim = imgs(i).n_sim(ind_y,:);
    fit_exp = imgs(i).fit_exp.fit(imgs(i).X,imgs(i).Y);
    fit_sim = imgs(i).fit_sim.fit(imgs(i).X,imgs(i).Y);
    tr(i).n_exp_fit = fit_exp(ind_y,:);
    tr(i).n_sim_fit = fit_sim(ind_y,:);

    tr(i).v_sim_x = imgs(i).v_sim_x(ind_y,:);
    tr(i).v_sim_y = imgs(i).v_sim_y(:,ind_x)';
    tr(i).v_exp_x = 100.*interp1([os(i).local.x]./10,[os(i).local.vExp],tr(i).x,'linear',0);

    % fit velocities for effective temperature
    tau = @(sig,T) getTauExp(sig,T);
    gamma = @(t,sig,T) t./tau(sig,T).^2./(1 + t.^2./tau(sig,T).^2);
    v = @(r,t,sig,T) r.*gamma(t,sig,T);
    

    r_cut = 0.24;
    ind_x = abs(tr(i).x) < r_cut;
    ind_y = abs(tr(i).y) < r_cut;

    fit_model = @(c,xdata) v(xdata,t(i),sig(1).sim_x,c);
    tr(i).T_x = lsqcurvefit(fit_model,data.Te,tr(i).x(ind_x),tr(i).v_sim_x(ind_x),0,1e3);
    tr(i).v_sim_x_fit = fit_model(tr(i).T_x,tr(i).x);

    fit_model = @(c,xdata) v(xdata,t(i),sig(1).sim_y,c);
    tr(i).T_y = lsqcurvefit(fit_model,data.Te,tr(i).y(ind_y),tr(i).v_sim_y(ind_y),0,1e3);
    tr(i).v_sim_y_fit = fit_model(tr(i).T_y,tr(i).y);

    fit_model = @(c,xdata) v(xdata,t(i),sig(1).exp_x,c);
    tr(i).T_exp = lsqcurvefit(fit_model,data.Te,tr(i).x(ind_x),tr(i).v_exp_x(ind_x),0,1e3);
    tr(i).v_exp_x_fit = fit_model(tr(i).T_exp,tr(i).x);

end

%% evaluate kinetic model

km = struct;
[~,km_sig,gam,Ti,Te,n2,n3] = kinetic_model(data.grids.time,data.sig0,data.Ti,data.Te,imgs(1).fit_sim.amp,data.config.eic_opt);
for i = 1:length(data.grids.time)
    km(i).sigx = km_sig(i);
    km(i).sigy = km_sig(i);
    km(i).sig = km_sig(i);
    km(i).Ti = Ti(i);
    km(i).Te = Te(i);
    km(i).n2 = n2(i);
    km(i).n3 = n3(i);
    km(i).vx = gam(i).*data.grids.x_vec;
    km(i).vy = gam(i).*data.grids.y_vec;
end

%% fit gaussain to all time points
for i = 1:length(data.grids.time)    
    x = data.grids.x_vec;
    y = data.grids.y_vec;
    if strcmp(data.config.eq_set,'ideal_2F')
        img = data.grids.vars(i).i_n;
    else
        img = data.grids.vars(i).n;
    end
    [gauss_fit(i)] = fitImgWithGaussian(x,y,img,0);
    disp(['2D Gaussian Fits: ' num2str(i) '/' num2str(length(data.grids.time))])
end

%% fit plasma sizes to extract effective temperatures
sig_fit = struct();

% fit x axis of simulation
tau = @(sig,T) getTauExp(sig,T);
fitmodel = @(c,data) c(1)*sqrt(1+data.^2/tau(c(1),c(2))^2);
p0 = [sig(1).sim_x os(1).Te];
xdata = t;
zdata = [sig.sim_x];
[c] = lsqcurvefit(fitmodel,p0,xdata,zdata,0,200);
sig_fit.sim_x = fitmodel(c,xdata);
sig_fit.sim_x_T = c(2);

% fit x axis of experiment
tau = @(sig,T) getTauExp(sig,T);
fitmodel = @(c,data) c(1)*sqrt(1+data.^2/tau(c(1),c(2))^2);
p0 = [sig(1).exp_x os(1).Te];
xdata = t;
zdata = [sig.exp_x];
[c] = lsqcurvefit(fitmodel,p0,xdata,zdata,0,200);
sig_fit.exp_x = fitmodel(c,xdata);
sig_fit.exp_x_T = c(2);

% fit y axis of simulation
tau = @(sig,T) getTauExp(sig,T);
fitmodel = @(c,data) c(1)*sqrt(1+data.^2/tau(c(1),c(2))^2);
p0 = [sig(1).sim_y os(1).Te];
xdata = t;
zdata = [sig.sim_y];
[c] = lsqcurvefit(fitmodel,p0,xdata,zdata,0,200);
sig_fit.sim_y = fitmodel(c,xdata);
sig_fit.sim_y_T = c(2);

% fit x axis of experiment
tau = @(sig,T) getTauExp(sig,T);
fitmodel = @(c,data) c(1)*sqrt(1+data.^2/tau(c(1),c(2))^2);
p0 = [sig(1).exp_y os(1).Te];
xdata = t;
zdata = [sig.exp_y];
[c] = lsqcurvefit(fitmodel,p0,xdata,zdata,0,200);
sig_fit.exp_y = fitmodel(c,xdata);
sig_fit.exp_y_T = c(2);

% fit 2D axis of experiment
tau = @(sig,T) getTauExp(sig,T);
fitmodel = @(c,data) c(1)*sqrt(1+data.^2/tau(c(1),c(2))^2);
p0 = [sig(1).exp_2D os(1).Te];
xdata = t;
zdata = [sig.exp_2D];
[c] = lsqcurvefit(fitmodel,p0,xdata,zdata,0,200);
sig_fit.exp_2D = fitmodel(c,xdata);
sig_fit.exp_2D_T = c(2);

% fit 2D axis of simulation
tau = @(sig,T) getTauExp(sig,T);
fitmodel = @(c,data) c(1)*sqrt(1+data.^2/tau(c(1),c(2))^2);
p0 = [sig(1).sim_2D os(1).Te];
xdata = t;
zdata = [sig.sim_2D];
[c] = lsqcurvefit(fitmodel,p0,xdata,zdata,0,200);
sig_fit.sim_2D = fitmodel(c,xdata);
sig_fit.sim_2D_T = c(2);

%% plot density images
% plot density with residuals
colvar = {'n_sim','n_exp','n_res'};
colstr = {'Sim','Exp','Sim - Exp'};
num = length(colvar);
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);

for k = 1:length(t)
    disp(['Plotting: ' num2str(k) '/' num2str(length(t))])
    
    % update axis information
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end
            
            cax = get_axis(fig,ax{i,j});
            xdata = imgs(k).x;
            ydata = imgs(k).y;
            zdata = imgs(k).(colvar{j});
            imagesc(xdata,ydata,zdata)
            colorbar
            cax.YDir = 'Normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            if i == size(ax,1), xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title(colstr{j},'FontWeight','normal')
            if j==1, cmax = max(zdata,[],'all'); end
            cax.CLim = [0 cmax*1.1];
            if j==3, cax.CLim = [-cmax cmax]/2; end

        end
    end
    str1 = ['t = ' num2str(t(k)*1e6,'%.3g') '\mus = ' num2str(t(k)/data.tau,'%.3g') '\tau_e_x_p'];
    an.String = str1;
    an.Position = [0.159500000000000,0.929069291338582,0.723000000000000,0.080100000000000];
    saveas(fig,[data.folder f 'imgs-' num2str(k) '.png']);
end
close(fig)

%% plot transect data

vars = {'v_sim_x','v_exp_x','v_sim_y','v_exp_x_fit'};
num = length(tr);
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
iter = 0;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end
        
        cax = get_axis(fig,ax{i,j});
        hold on
        l = get_line_specs(length(vars));
        for k = 1:length(vars)
            xdata = tr(iter).x;
            if contains(vars{k},'_y')
                xdata = tr(iter).y;
            end
            ydata = tr(iter).(vars{k});
            plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(k).col,'MarkerFaceColor',l(k).col,'MarkerEdgeColor',l(k).col)
        end
        xlabel('x (cm)')
        ylabel('n (cm^-^3)')
        titleStr = ['t = ' num2str(t(iter)*1e6,'%.3g') ' \mus'];
        title(titleStr)

    end
end
lgd = legend(vars);
lgd.Position = [0.899896476035170,0.198423988232002,0.080329755350856,0.236471041927977];
saveas(fig,[data.folder f 'velocity-fits.png']);
close(fig)

%% plot density and velocity transects together

% select which time points to plot
t_plot = linspace(0,max(t),4);
ind_t = zeros(size(t_plot));
for i = 1:length(t_plot)
    [~,ind_t(i)] = min(abs(t - t_plot(i)));
end
t_plot = t(ind_t);

% create figure
vars{1} = {'n_exp','n_sim'};
vars{2} = {'v_exp_x','v_sim_x'};
row_str = {'n (10^8 cm^-^3)','v (cm/s)'};
lgd_str = {'data','sim'};
v_cut = [0.45 0.45 0.45 0.45];

num = 8;
[fig,ax,an] = open_subplot(num);
iter = 0;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end
        
        cax = get_axis(fig,ax{i,j});

        hold off
        l = get_line_specs(length(vars{i}));
        for k = 1:length(vars{i})
            xdata = [tr(ind_t(j)).x];
            ydata = [tr(ind_t(j)).(vars{i}{k})];
            if i == 2
                ind_x = abs(xdata) < v_cut(j);
                xdata = xdata(ind_x);
                ydata = ydata(ind_x);
            end
            if i == 1, ydata = ydata./1e8; end
            plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(k).col,'MarkerFaceColor',l(k).col,'MarkerEdgeColor',l(k).col)
            hold on
        end

        cax.FontSize = 11;
        if i == size(ax,1), xlabel('x (cm)'), end
        if j == 1, ylabel(row_str{i}), end
        if i == 1, title([num2str(t_plot(j)*1e6) ' \mus'],'FontWeight','normal'); end
        xlim([min(tr(ind_t(j)).x) max(tr(ind_t(j)).x)])
    end
end
an.String = [];
lgd = legend(lgd_str);
lgd.Position = [0.825186396355424,0.127480879107887,0.075254105375184,0.091235634102218];
saveas(fig,[data.folder f 'n-v-transects.png']);
close(fig)



%% create size evolution figure for paper

ind_t = data.grids.time < max(t);
Te = zeros(size(data.grids.time(ind_t)));
Ti = zeros(size(data.grids.time(ind_t)));
for i = 1:length(data.grids.time(ind_t))
    [~,ind_x] = min(abs(data.grids.x_vec));
    [~,ind_y] = min(abs(data.grids.y_vec));
    Te(i) = data.grids.vars(i).e_temp(ind_y,ind_x);
    Ti(i) = data.grids.vars(i).i_temp(ind_y,ind_x);
end

Ti_exp = zeros(size(os));
for i = 1:length(os)
    x = [os(i).local.x]./10;
    ind_x = abs(x) < sig(1).exp_x/2;
    Ti_exp(i) = mean([os(i).local(ind_x).Ti]);
end

num = 4;
[fig,ax,an] = open_subplot(num);

cax = get_axis(fig,ax{1,1});
l = get_line_specs(3);
lp = l(1);
plot(t*1e6,[sig.exp_x],'o','LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
hold on
lp = l(2);
plot(data.grids.time(ind_t)*1e6,[gauss_fit(ind_t).sigx],'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
hold on
lp = l(3);
plot(data.grids.time(ind_t)*1e6,[km(ind_t).sig],'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
xlabel('t (\mus)')
ylabel('\sigma_x (cm)')

cax = get_axis(fig,ax{1,2});
l = get_line_specs(3);
lp = l(1);
plot(t*1e6,[sig.exp_y],'o','LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
hold on
lp = l(2);
plot(data.grids.time(ind_t)*1e6,[gauss_fit(ind_t).sigy],'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
hold on
lp = l(3);
plot(data.grids.time(ind_t)*1e6,[km(ind_t).sig],'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
xlabel('t (\mus)')
ylabel('\sigma_y (cm)')

lgd = legend({'exp','sim','Vlasov'});
lgd.Position = [0.765424889057839,0.597624081941271,0.125958703274572,0.085159012420438];

cax = get_axis(fig,ax{2,1});
l = get_line_specs(3);
lp = l(2);
plot(data.grids.time(ind_t)*1e6,Te,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
hold on
lp = l(3);
plot(data.grids.time(ind_t)*1e6,[km(ind_t).Te],'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
xlabel('t (\mus)')
ylabel('T_e (K)')


cax = get_axis(fig,ax{2,2});
l = get_line_specs(3);
hold on
lp = l(1);
plot(t*1e6,Ti_exp,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
lp = l(2);
plot(data.grids.time(ind_t)*1e6,Ti,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
lp = l(3);
plot(data.grids.time(ind_t)*1e6,[km(ind_t).Ti],'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
xlabel('t (\mus)')
ylabel('T_i (K)')

saveas(fig,[data.folder f 'size-temp-evol.png']);
close(fig)

%% output results
s.imgs = imgs;
s.tr = tr;
s.sig = sig;
s.sig_fit = sig_fit;
s.t = t;
end