function [] = iawAnalysis(data,flags)
time_format = '%.3g';
%% Initialize IAW Data Struct
iaw = struct();
num_time_points = length(data.grids.time);
iaw(num_time_points).t = [];
for i = 1:length(iaw)
    iaw(i).t = data.grids.time(i); % cgs units
end

%% Extract 1D Density Transects from 2D Grids
for i = 1:length(iaw)
    y_cut_off_distance = data.grids.gauss_fits(i).sigy/2;
    y_points_to_average = abs(data.grids.y_vec) < y_cut_off_distance;

    iaw(i).x = data.grids.x_vec;
    iaw(i).n = mean(data.grids.vars(i).n(y_points_to_average,:),1);
    iaw(i).Ti = mean(data.grids.vars(i).i_temp(y_points_to_average,:),1);
    iaw(i).Te = mean(data.grids.vars(i).e_temp(y_points_to_average,:),1);
end

%% Fit 1D Gaussian to Density Distribution

fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
frames = cell(1,length(iaw));
grid minor
for i = 1:length(iaw)
    iaw(i).n_gfit = fitGaussianNoOffset1D(iaw(i).x,iaw(i).n);
    iaw(i).gsig = iaw(i).n_gfit.sig;
    cla
    hold on
    plot(iaw(i).x,iaw(i).n,'.','MarkerSize',8)
    plot(iaw(i).x,iaw(i).n_gfit.guess,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).n_gfit.fit,'LineWidth',1.5)
    legend({'data','guess','fit'})
    xlabel('x (cm)')
    ylabel('n (10^8 cm^-^3)')
    title(['t = ' num2str(iaw(i).t*1e6,time_format) '\mus'],'FontWeight','normal')
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'iaw-gauss-fits'],frames)

%% Fit 1D Gaussian with Slope to Density Distribution

fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
frames = cell(1,length(iaw));
grid minor
for i = 1:length(iaw)
    iaw(i).n_gfit_slope = fitGaussianNoOffsetWithSlope1D(iaw(i).x,iaw(i).n);
    
    cla
    hold on
    plot(iaw(i).x,iaw(i).n,'.','MarkerSize',8)
    plot(iaw(i).x,iaw(i).n_gfit_slope.guess,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).n_gfit_slope.fit,'LineWidth',1.5)
    legend({'data','guess','fit'})
    xlabel('x (cm)')
    ylabel('n (10^8 cm^-^3)')
    title(['t = ' num2str(iaw(i).t*1e6,time_format) '\mus'],'FontWeight','normal')
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'iaw-gauss-fits-slope'],frames)

%% Fit 1D Gaussian with Slope Offset to Density Distribution

fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
frames = cell(1,length(iaw));
grid minor
for i = 1:length(iaw)
    iaw(i).n_gfit_slope_offset = fitGaussianNoOffsetWithSlopeOffset1D(iaw(i).x,iaw(i).n);
    
    
    cla
    hold on
    plot(iaw(i).x,iaw(i).n,'.','MarkerSize',8)
    plot(iaw(i).x,iaw(i).n_gfit_slope_offset.guess,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).n_gfit_slope_offset.fit,'LineWidth',1.5)
    legend({'data','guess','fit'})
    xlabel('x (cm)')
    ylabel('n (10^8 cm^-^3)')
    title(['t = ' num2str(iaw(i).t*1e6,time_format) '\mus'],'FontWeight','normal')
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'iaw-gauss-fits-slope-offset'],frames)

%% Compute Density Perturbation from Gauss Fits

fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
frames = cell(1,length(iaw));
grid minor
for i = 1:length(iaw)
    iaw(i).dn = iaw(i).n - iaw(i).n_gfit.fit;
    iaw(i).dntild = iaw(i).dn./iaw(i).n_gfit.amp;

    iaw(i).dn_slope = iaw(i).n - iaw(i).n_gfit_slope.fit;
    iaw(i).dntild_slope = iaw(i).dn_slope./iaw(i).n_gfit_slope.amp;

    iaw(i).dn_slope_offset = iaw(i).n - iaw(i).n_gfit_slope_offset.fit;
    iaw(i).dntild_slope_offset = iaw(i).dn_slope_offset./iaw(i).n_gfit_slope_offset.amp;

    cla
    hold on
    plot(iaw(i).x,iaw(i).dntild,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_slope,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_slope_offset,'LineWidth',1.5)
    legend({'dntild','dntild_s','dntild_so'})
    xlabel('x (cm)')
    ylabel('\deltan(x,t)/n(0,t)')
    title(['t = ' num2str(iaw(i).t*1e6,time_format) '\mus'],'FontWeight','normal')
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'iaw-dntild'],frames)

%% Fit Density Perturbations with IAW Model

fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
frames = cell(1,length(iaw));
grid minor
for i = 1:length(iaw)
    iaw(i).tau = getTauExp(iaw(i).n_gfit.sig,data.Te+data.Ti);
    iaw(i).time_fac = sqrt(1+iaw(i).t^2/iaw(i).tau^2);
    iaw(i).k_vlasov = 2*pi/(data.settings.n_iaw_sig*iaw(i).time_fac);
    iaw(i).dntild_fit = fitIAW1D(iaw(i).x,iaw(i).dntild,iaw(i).k_vlasov);

    cla
    hold on
    plot(iaw(i).x,iaw(i).dntild,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_fit.guess,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_fit.fit,'LineWidth',1.5)
    legend({'data','guess','fit'})
    xlabel('x (cm)')
    ylabel('\deltan(x,t)/n(0,t)')
    title(['t = ' num2str(iaw(i).t*1e6,time_format) '\mus'],'FontWeight','normal')
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'iaw-dntild-fits'],frames)

%% Fit Density Perturbations with IAW Model

fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
frames = cell(1,length(iaw));
grid minor
for i = 1:length(iaw)
    iaw(i).tau = getTauExp(iaw(i).n_gfit.sig,data.Te+data.Ti);
    iaw(i).time_fac = sqrt(1+iaw(i).t^2/iaw(i).tau^2);
    iaw(i).k_vlasov = 2*pi/(data.settings.n_iaw_sig*iaw(i).time_fac);
    iaw(i).dntild_slope_fit = fitIAW1D(iaw(i).x,iaw(i).dntild_slope,iaw(i).k_vlasov);

    cla
    hold on
    plot(iaw(i).x,iaw(i).dntild,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_fit.guess,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_fit.fit,'LineWidth',1.5)
    legend({'data','guess','fit'})
    xlabel('x (cm)')
    ylabel('\deltan(x,t)/n(0,t)')
    title(['t = ' num2str(iaw(i).t*1e6,time_format) '\mus'],'FontWeight','normal')
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'iaw-dntild-slope-fits'],frames)

%% Fit Density Perturbations with IAW Model

fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
frames = cell(1,length(iaw));
grid minor
for i = 1:length(iaw)
    iaw(i).tau = getTauExp(iaw(1).n_gfit.sig,data.Te+data.Ti);
    iaw(i).time_fac = sqrt(1+iaw(i).t^2/iaw(i).tau^2);
    iaw(i).k_vlasov = 2*pi/(data.settings.n_iaw_sig*iaw(i).time_fac);
    iaw(i).dntild_slope_offset_fit = fitIAW1D(iaw(i).x,iaw(i).dntild_slope_offset,iaw(i).k_vlasov);

    cla
    hold on
    plot(iaw(i).x,iaw(i).dntild,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_fit.guess,'LineWidth',1.5)
    plot(iaw(i).x,iaw(i).dntild_fit.fit,'LineWidth',1.5)
    legend({'data','guess','fit'})
    xlabel('x (cm)')
    ylabel('\deltan(x,t)/n(0,t)')
    title(['t = ' num2str(iaw(i).t*1e6,time_format) '\mus'],'FontWeight','normal')
    frames{i} = getframe(fig);
end
close(fig)
write_video([data.folder filesep 'iaw-dntild-slope-offset-fits'],frames)

%% Gather Fit Parameters
for i = 1:length(iaw)
    iaw(i).A = iaw(i).dntild_fit.amp;
    iaw(i).A_slope = iaw(i).dntild_slope_fit.amp;
    iaw(i).A_slope_offset = iaw(i).dntild_slope_offset_fit.amp;
    iaw(i).lam = iaw(i).dntild_fit.lam;
    iaw(i).lam_slope = iaw(i).dntild_slope_fit.lam;
    iaw(i).lam_slope_offset = iaw(i).dntild_slope_offset_fit.lam;
    iaw(i).k = iaw(i).dntild_fit.k;
    iaw(i).k_slope = iaw(i).dntild_slope_fit.k;
    iaw(i).k_slope_offset = iaw(i).dntild_slope_offset_fit.k;
    iaw(i).sig = iaw(i).dntild_fit.g_sig;
    iaw(i).sig_slope = iaw(i).dntild_slope_fit.g_sig;
    iaw(i).sig_slope_offset = iaw(i).dntild_slope_offset_fit.g_sig;
end

%% Plot Wavelength Fit Parameters
fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
cla
hold on
plot([iaw.t],[iaw.lam],'.','MarkerSize',8)
plot([iaw.t],[iaw.lam_slope],'.','MarkerSize',8)
plot([iaw.t],[iaw.lam_slope_offset],'.','MarkerSize',8)
legend({'G','G_s','G_s_o'})
xlabel('time (\mus)')
ylabel('\lambda (cm)')
grid minor
exportgraphics(fig,[data.folder filesep 'iaw-lam.png'],'Resolution',300)
close(fig)

%% Plot Amplitude Fit Parameters
fig = figure;
fig.Position = [198.6000  213.8000  918.4000  438.8000];
cla
hold on
plot([iaw.t],[iaw.A],'.','MarkerSize',8)
plot([iaw.t],[iaw.A_slope],'.','MarkerSize',8)
plot([iaw.t],[iaw.A_slope_offset],'.','MarkerSize',8)
legend({'G','G_s','G_s_o'})
xlabel('time (\mus)')
ylabel('A(t)')
grid minor
exportgraphics(fig,[data.folder filesep 'iaw-amp.png'],'Resolution',300)
close(fig)

%% Compute Dispersion Relations

for i = 1:length(iaw)
    iaw(i).Te_global = mean(iaw(i).Te,'all');
    iaw(i).Ti_global = trapz(iaw(i).x,iaw(i).n.*iaw(i).Ti)/trapz(iaw(i).x,iaw(i).n);
    iaw(i).cs = getSoundSpeed(data.settings.m_i,iaw(i).Ti_global,iaw(i).Te_global,data.adiabatic_index_e,data.adiabatic_index_i);
    iaw(i).w = iaw(i).k*iaw(i).cs;
    iaw(i).w_slope = iaw(i).k_slope*iaw(i).cs;
    iaw(i).w_slope_offset = iaw(i).k_slope_offset*iaw(i).cs;
end

fit_model = @(p,t) p(2)./(1+t.^2./p(1).^2);
xdata = [iaw.t];
ydata = [iaw.k_slope_offset].*sqrt([iaw.Te_global]);
p0 = [iaw(1).tau ydata(1)];
p = lsqcurvefit(fit_model,p0,xdata,ydata);
% amp_fit = p(1);
tau_fit = p(1);

fig = figure;
hold on
plot(xdata.*1e6,ydata,'.','MarkerSize',12)
plot(xdata.*1e6,fit_model(p,xdata),'LineWidth',1.5)
xlabel('t (\mus)')
ylabel('$k\sqrt{T_e}$','Interpreter','latex')
title(['Fit: \tau = ' num2str(tau_fit*1e6,time_format) '\mus'],'FontWeight','normal')
exportgraphics(fig,[data.folder filesep 'iaw-fit-kRootTe.png'],'Resolution',300)
close(fig)

%% Fit Amplitude Evolution
disp_fit = fitIAWDispersion([iaw.t],[iaw.A_slope_offset],iaw(1).w_slope_offset,tau_fit);

fig = figure;
hold on
plot([iaw.t],[iaw.A_slope_offset],'.','MarkerSize',12)
plot([iaw.t],disp_fit.fit,'-','LineWidth',1.5)
legend('data','fit')
xlabel('t (\mus)')
ylabel('A(t)')
title('Data Analyzed with Gaussian + Slope Offset','FontWeight','normal')
exportgraphics(fig,[data.folder filesep 'iaw-fit-to-A.png'],'Resolution',300)
close(fig)

%% Save IAW Struct
save([data.folder filesep 'iaw.mat'],'iaw','disp_fit');

end