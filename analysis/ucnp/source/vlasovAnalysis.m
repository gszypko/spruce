function [] = vlasovAnalysis(data,flags)
%% Ensure 2D Gaussian Fits are Done

if ~flags.doGaussianFits2D, doGaussianFits2D(data,flags); end
load([data.folder filesep 'gauss-fits.mat'],'gauss_fits');

%% Compile MHD Data and Vlasov Theory for Comparison
vlasov = struct;

vlasov.sim = struct;
vlasov.time = data.grids.time;
for i = 1:length(data.grids.time)
    
    vlasov.sim(i).sigx = gauss_fits(i).sigx;
    vlasov.sim(i).sigy = gauss_fits(i).sigy;
    vlasov.sim(i).sig = sqrt(gauss_fits(i).sigx*gauss_fits(i).sigy);
    vlasov.sim(i).n2 = gauss_fits(i).amp;
    vlasov.sim(i).n3 = gauss_fits(i).amp;
    
    [~,indx] = min(abs(data.grids.x_vec));
    [~,indy] = min(abs(data.grids.y_vec));
        
    vlasov.sim(i).Ti = data.grids.vars(i).(data.i.T)(indy,indx);
    vlasov.sim(i).Te = data.grids.vars(i).(data.e.T)(indy,indx);
    vlasov.sim(i).T = vlasov.sim(i).Ti + vlasov.sim(i).Te;
    
    vlasov.sim(i).xt.r = data.grids.pos_x(indy,:);
    vlasov.sim(i).yt.r = data.grids.pos_y(:,indx)';
    vlasov.sim(i).xt.n = data.grids.vars(i).(data.i.n)(indy,:);
    vlasov.sim(i).yt.n = data.grids.vars(i).(data.i.n)(:,indx)';
    vlasov.sim(i).xt.v = data.grids.vars(i).(data.i.v_x)(indy,:);
    vlasov.sim(i).yt.v = data.grids.vars(i).(data.i.v_y)(:,indx)';
    vlasov.sim(i).xt.Te = data.grids.vars(i).(data.e.T)(indy,:);
    vlasov.sim(i).yt.Te = data.grids.vars(i).(data.e.T)(:,indx)';
    vlasov.sim(i).xt.Ti = data.grids.vars(i).(data.i.T)(indy,:);
    vlasov.sim(i).yt.Ti = data.grids.vars(i).(data.i.T)(:,indx)';
    vlasov.sim(i).xt.T = vlasov.sim(i).yt.Te + vlasov.sim(i).yt.Ti;
    vlasov.sim(i).yt.T = vlasov.sim(i).yt.Te + vlasov.sim(i).yt.Ti;
    
end
vlasov.tau = getTauExp(vlasov.sim(1).sig,vlasov.sim(1).T);

vlasov.theory = struct;
[~,sig,gam,Ti,Te,n2,n3] = kinetic_model(vlasov.time,vlasov.sim(1).sig,vlasov.sim(1).Ti,vlasov.sim(1).Te,vlasov.sim(1).n2,data.config.eic_opt);
for i = 1:length([vlasov.time])
    vlasov.theory(i).sigx = sig(i);
    vlasov.theory(i).sigy = sig(i);
    vlasov.theory(i).sig = sig(i);
    vlasov.theory(i).Ti = Ti(i);
    vlasov.theory(i).Te = Te(i);
    vlasov.theory(i).T = Ti(i) + Te(i);
    vlasov.theory(i).n2 = n2(i);
    vlasov.theory(i).n3 = n3(i);
    vlasov.theory(i).gam = gam(i);
    vlasov.theory(i).v_x = gam(i).*vlasov.sim(i).xt.r;
    vlasov.theory(i).v_y = gam(i).*vlasov.sim(i).yt.r;
end

%% Fit Size Evolution for Effective Temperature

% define functions for velocity along a given axis
tau = @(sig,T) getTauExp(sig,T);
sig = @(t,sig0,T) sig0.*sqrt(1+t.^2./tau(sig0,T).^2);
gam = @(t,sig,T) t./tau(sig,T).^2/(1+t.^2./tau(sig,T).^2);
v = @(r,t,sig,T) r.*gam(t,sig,T);

% fit for temperature on x axis
fit_model_x = @(c,data) sig(data,vlasov.sim(1).sigx,c);
xdata = [vlasov.time];
ydata = [vlasov.sim.sigx];
T_x = lsqcurvefit(fit_model_x,vlasov.sim(1).T,xdata,ydata,0,1e3);
vlasov.T_x = T_x;

% fit for temperature on y axis
fit_model_y = @(c,data) sig(data,vlasov.sim(1).sigy,c);
xdata = [vlasov.time];
ydata = [vlasov.sim.sigy];
T_y = lsqcurvefit(fit_model_y,vlasov.sim(1).T,xdata,ydata,0,1e3);
vlasov.T_y = T_y;

% fit for rms temperature
fit_model_2D = @(c,data) sig(data,vlasov.sim(1).sig,c);
xdata = [vlasov.time];
ydata = [vlasov.sim.sig];
T_2D = lsqcurvefit(fit_model_2D,vlasov.sim(1).T,xdata,ydata,0,1e3);
vlasov.T_2D = T_2D;

% store fits in data structure
sigx_fit = fit_model_x(vlasov.T_x,vlasov.time);
sigy_fit = fit_model_y(vlasov.T_y,vlasov.time);
sig_fit = fit_model_2D(vlasov.T_2D,vlasov.time);
for i = 1:length(vlasov.time)
    vlasov.sim(i).sigx_fit = sigx_fit(i);
    vlasov.sim(i).sigy_fit = sigy_fit(i);
    vlasov.sim(i).sig_fit = sig_fit(i);
end

% plot fits to plasma size
p.x.data = [vlasov.sim.sigx];
p.x.fit = [vlasov.sim.sigx_fit];
p.y.data = [vlasov.sim.sigy];
p.y.fit = [vlasov.sim.sigy_fit];
p.gm.data = [vlasov.sim.sig];
p.gm.fit = [vlasov.sim.sig_fit];
fields = fieldnames(p);

colvar = {'x','y','gm'};
colstr = {'\sigma_x (cm)','\sigma_y (cm)','\sigma (cm)'};
vars = {'data','fit'};
title_var = {'T_x','T_y','T_2D'};

num = length(fields);
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        cax = get_axis(fig,ax{i,j});
        hold on

        l = get_line_specs(length(vars));
        for k = 1:length(vars)
            xdata = [vlasov.time].*1e6;
            ydata = [p.(colvar{iter}).(vars{k})];
            lp = l(k);
            plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
        end
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 12;
        grid on
        grid minor

        if i == row, xlabel('t (\mus)'), end
        ylabel(colstr{iter})
        title(['T_e = ' num2str(vlasov.(title_var{iter})) ' K'],'FontWeight','normal')
    end
end
lgd = legend(vars);
lgd.Position = [0.814763192769338,0.212150191817874,0.084684301063469,0.152255642324462];
saveas(fig,[data.folder filesep 'sig-fits.png']);
close(fig)

%% Plot Summary of Expansion for MHD Sims Against Vlasov Theory

q = {'sim','theory'};
qstr{1} = 'Simulation';
if data.config.eic_opt, qstr{2} = 'Vlasov';
else, qstr{2} = 'Vlasov'; end
colvar = {'sigx','sigy','sig','n2','Ti','Te'};
colstr = {'\sigma_x (cm)','\sigma_y (cm)','\sigma (cm)','n_0 (10^8 cm^-^3)','T_i (K)','T_e (K)'};
num = length(colvar);
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        cax = get_axis(fig,ax{i,j});
        hold on

        lgdstr = qstr;
        l = get_line_specs(length(lgdstr));
        
        for k = 1:length(q)
            xdata = [vlasov.time]./vlasov.tau;
            ydata = [vlasov.(q{k}).(colvar{iter})];
            lp = l(k);
            plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
        end
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 12;
        grid on
        grid minor

        if i == row, xlabel('t / \tau'), end
        ylabel(colstr{iter})
    end
end

% add legend
lgd = legend(lgdstr);
lgd.Position = [0.8954    0.4610    0.0911    0.0787];

saveas(fig,[data.folder filesep 'vlasov-evol.png']);
close(fig)

%% Plot Central Transects for Density, Velocity, and Temperature

% generate figure
frames = cell(length([vlasov.time]));
num = 6;
[fig,ax,an,row,col] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.191979522184300,0.890899999999999,0.687720477815701,0.080100000000000];

col_var = {'n','v','Te'};
row_var_y = {'xt','yt'};
col_str = {'n (cm^-3)','v (cm/s)','T_e (K)'};
row_str = {'x (cm)','y (cm)'};

for k = 1:length(vlasov.time)
    disp(['Plotting Transects: ' num2str(k) '/' num2str(length(vlasov.time))])
    hold([ax{:}],'off')

    iter = 0;
    for i = 1:row
        for j = 1:col
            if iter > num - 1, break, end
            iter = iter + 1;
            cax = get_axis(fig,ax{i,j});
            hold off
            
            l = get_line_specs(2);
            for indq = 1:length(q)
                xdata = [vlasov.sim(k).(row_var_y{i}).r];
                ydata = [vlasov.sim(k).(row_var_y{i}).(col_var{j})];
                lp = l(1);
                plot(xdata,ydata,lp.style,'LineWidth',2,'MarkerSize',4,'Color',lp.col,'MarkerFaceColor',lp.col,'MarkerEdgeColor',lp.col)
            end
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            grid on
            grid minor
            xlim([-.4 .4])
            xlabel(row_str{i})
            ylabel(col_str{j})

        end
    end


    an.String = ['t = ' num2str(vlasov.time(k)*1e6,'%.2g') ' \mus = ' num2str(vlasov.time(k)/vlasov.tau,'%.2g') ' \tau_{exp}'];
    frames{k} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'transects'],frames);

save([data.folder filesep 'vlasov.mat'],'vlasov','-mat');

end