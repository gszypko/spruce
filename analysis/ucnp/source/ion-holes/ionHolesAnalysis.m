function [] = ionHolesAnalysis(data,flags)
%% Check for Required Grids
required_grids = {'n','v_x','v_y','i_temp','e_temp'};
vars = fieldnames(data.grids.vars);
for i = 1:length(required_grids)
    if ~max(strcmp(vars,required_grids{i}))
        error(['Expected variable <' required_grids{i} '> is missing.'])
    end
end

%% Load Experimental Data - If Applicable
% check if <os.mat> is within parent folder
f = filesep;
parent_folder = extractBefore(data.folder,[f 'set']);
files = dir(parent_folder);
file_ind = strcmp({files.name},'os.mat');
exp_file_found = max(file_ind);

% if <os.mat> is found, load it and adjust delay times to include half of exposure time + convert to cgs units
if exp_file_found
    load([parent_folder f 'os.mat'],'os');
    for i = 1:length(os)
        os(i).delays = (os(i).delays + os(i).tE/2)*1e-9;
    end
end
os = os([os.delays] < 22e-6);

% find indices of simulation time points that are closest to experimental time points
% ind_t(i) holds the index of data.grids.vars corresponding to experimental time point os(i).delays
closest_time_point = zeros(size(os)); 
for i = 1:length([os.delays])
    [~,closest_time_point(i)] = min(abs(os(i).delays - data.grids.time));
end

% select which time points to plot
t = [os.delays];
t_plot = [os.delays];
ind_t_plot = zeros(size(t_plot));
for i = 1:length(t_plot)
    [~,ind_t_plot(i)] = min(abs(t - t_plot(i)));
end
t_plot = t(ind_t_plot);

%% Create Structure of Information

% set aside relevant grid info for simulation
s = struct();
s.folder = data.folder;
s.sim_grids(length(data.grids.time)).t = [];
for i = 1:length(data.grids.time)
    ind_x = data.grids.x_vec > -flags.plot_window(1) & data.grids.x_vec < flags.plot_window(1);
    ind_y = data.grids.y_vec > -flags.plot_window(2) & data.grids.y_vec < flags.plot_window(2);

    s.sim_grids(i).t = data.grids.time(i);
    s.sim_grids(i).x = data.grids.x_vec(ind_x);
    s.sim_grids(i).y = data.grids.y_vec(ind_y);
    s.sim_grids(i).n = data.grids.vars(i).n(ind_y,ind_x);
    s.sim_grids(i).Ti = data.grids.vars(i).i_temp(ind_y,ind_x);
    s.sim_grids(i).Te = data.grids.vars(i).e_temp(ind_y,ind_x);
    s.sim_grids(i).v_x = data.grids.vars(i).v_x(ind_y,ind_x);
    s.sim_grids(i).v_y = data.grids.vars(i).v_y(ind_y,ind_x);
end

% set aside relevant grid info for experiment
if exp_file_found
    init_exp_n = os(1).imgs.density*1e8;
    init_exp_n_sg = sgfilt2D(init_exp_n,7,7,1,1,true);
    max_exp_n_sg = max(init_exp_n_sg,[],'all');
    s.exp_grids(length(os)).t = [];
    [~,~,vlasov_gam,vlasov_Ti,vlasov_Te,~,~] = kinetic_model([os.delays],data.sig0,data.Ti,data.Te,max_exp_n_sg,true);
    for i = 1:length(os)
        ind_x = os(i).imgs.xRelInMM/10 > -flags.plot_window(1) & os(i).imgs.xRelInMM/10 < flags.plot_window(1);
        ind_y = os(i).imgs.yRelInMM/10 > -flags.plot_window(2) & os(i).imgs.yRelInMM/10 < flags.plot_window(2);

        s.exp_grids(i).t = os(i).delays;
        s.exp_grids(i).x = os(i).imgs.xRelInMM(ind_x)./10;
        s.exp_grids(i).y = os(i).imgs.yRelInMM(ind_y)./10;
        s.exp_grids(i).n = os(i).imgs.density(ind_y,ind_x)*1e8;
        ones_mat = ones(size(s.exp_grids(i).n));
        s.exp_grids(i).Ti = vlasov_Ti(i).*ones_mat;
        s.exp_grids(i).Te = vlasov_Te(i).*ones_mat;
        s.exp_grids(i).v_x = vlasov_gam(i).*s.exp_grids(i).x.*ones_mat;
        s.exp_grids(i).v_y = vlasov_gam(i).*s.exp_grids(i).y'.*ones_mat;
    end
end


%% Determine Ion Hole Orientation
[~,ind_t] = min(data.grids.time);
x = data.grids.x_vec;
y = data.grids.y_vec;
n = data.grids.vars(ind_t).n;
[gauss_fit,Rx,Ry] = extractHoleOrientation(x,y,n,flags.figvis,flags.hole_orientation,data.folder);

s.adiabatic_index_i = data.adiabatic_index_i;
s.adiabatic_index_e = data.adiabatic_index_e;
s.m_i = data.settings.m_i;
s.Te = data.Te;
s.Ti = data.Ti;
s.sig_x0 = gauss_fit.sigx_g;
s.sig_y0 = gauss_fit.sigy_g;
s.tau_hole = getTauExp(gauss_fit.sigx_g,data.Te+data.Ti);

clearvars -except s Rx Ry gauss_fit exp_file_found t_plot closest_time_point ind_t_plot

%% Project Data on Grid Aligned with Ion Hole Propagation Axis

fields = fieldnames(s);
grid_ind = endsWith(fields,'_grids');
fields = fields(grid_ind);
for i = 1:length(fields)
    name = extractBefore(fields{i},'_grids');
    name_full = [name '_rot'];
    name_1D = [name '_1D'];
    name_int = [name '_int'];
    s2 = s.(fields{i});

    for j = 1:length(s2)
        x = s2(j).x;
        y = s2(j).y;
        [X,Y] = meshgrid(x,y);

        s.(name_full)(j).t = s2(j).t;

        x_rot = x - gauss_fit.x0;
        y_rot = y - gauss_fit.y0;
        [X_rot,Y_rot] = meshgrid(x_rot,y_rot);

        s.(name_full)(j).x = x;
        s.(name_full)(j).y = y;
        s.(name_full)(j).n = interp2(X,Y,s2(j).n,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
        s.(name_full)(j).Ti = interp2(X,Y,s2(j).Ti,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
        s.(name_full)(j).Te = interp2(X,Y,s2(j).Te,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
        v_x_rot = interp2(X,Y,s2(j).v_x,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
        v_y_rot = interp2(X,Y,s2(j).v_y,Rx(X_rot,Y_rot,-gauss_fit.theta),Ry(X_rot,Y_rot,-gauss_fit.theta));
        s.(name_full)(j).v = Rx(v_x_rot,v_y_rot,gauss_fit.theta);
        s.(name_full)(j).pi = s.m_i.*s.(name_full)(j).n.*s.(name_full)(j).v;
    end

    s.(name_1D) = struct();
    s.(name_1D)(length(s.(name_full))).t = [];
    for k = 1:length(s.(name_full))
        s.(name_1D)(k).t = s.(name_full)(k).t;
        time_fac = sqrt(1+s.(name_1D)(k).t.^2/s.tau_hole^2);
        ind_x = abs(s.(name_full)(k).x) < 3*s.sig_x0;
        ind_y = abs(s.(name_full)(k).y) < gauss_fit.sigy_g*time_fac/2;
        fields3 = fieldnames(s.(name_full));
        s.(name_1D)(k).r = s.(name_full)(k).x(ind_x);
        for l = 4:length(fields3)
            s.(name_1D)(k).(fields3{l}) = mean(s.(name_full)(k).(fields3{l})(ind_y,ind_x),1,'omitnan'); 
        end
    end
    s.(name_1D)(1).t = 0;

    % create anonymous functions to interpolate 1D variables in time and space
    s.(name_int) = struct();
    [R,T] = meshgrid(s.(name_1D)(1).r,[s.(name_1D).t]);
    n = reshape([s.(name_1D).n]',[],length(s.(name_1D)))';
    Te = reshape([s.(name_1D).Te]',[],length(s.(name_1D)))';
    Ti = reshape([s.(name_1D).Ti]',[],length(s.(name_1D)))';
    v = reshape([s.(name_1D).v]',[],length(s.(name_1D)))';
    s.(name_int).n = @(t,r) interp2(R,T,n,r,t);
    s.(name_int).v = @(t,r) interp2(R,T,v,r,t);
    s.(name_int).Te = @(t,r) interp2(R,T,Te,r,t);
    s.(name_int).Ti = @(t,r) interp2(R,T,Ti,r,t);
    s.(name_int).hole = @(t) ionHoleProp(t,s.(name_int).v,s.(name_int).Ti,s.(name_int).Te,s.m_i,s.adiabatic_index_e,s.adiabatic_index_i);

end

%% Fit 1D Data to Extract Ion Hole Location
fields = fieldnames(s);
ind_1D = endsWith(fields,'_1D');
ind_int = endsWith(fields,'_int');
fields_1D = fields(ind_1D);
fields_int = fields(ind_int);
fields = extractBefore(fields_int,'_int');
for i = 1:length(fields_1D)
    name_full = [fields{i} '_hole'];
    for k = 1:length(s.(fields_1D{i}))
        disp(['Fitting for ion hole location: ' num2str(k) '/' num2str(length(s.(fields_1D{i})))])
        s.(name_full)(k) = extractHoleLocation(s.(fields_1D{i})(k).t,s.(fields_1D{i})(k).r,s.(fields_1D{i})(k).n,s.(fields_1D{i})(k).v,s.sig_x0,s.tau_hole,s.(fields_int{i}).hole(s.(fields_1D{i})(k).t));
    end
end


%% Plot Transects with Hole Locations
fields = fieldnames(s);
ind_hole = endsWith(fields,'_hole');
ind_tau = ~contains(fields,'tau');
fields_hole = fields(ind_hole&ind_tau);

for l = 1:length(fields_hole)
    num = 4;
    [fig,ax,an] = open_subplot(num,'Visible','on');
    an.Position = [0.1567    0.9611    0.7230    0.0339];
    frames = cell(size(s.(fields_1D{l})));
    q = {'n','v','pi','Te'};
    qstr = {'n','v','\pi','T_e'};
    
    % update axis information
    for k = 1:length(s.(fields_hole{l}))
        disp(['Plotting transect fits: ' num2str(k) '/' num2str(length([s.(fields_1D{l}).t]))])
        iter = 0;
        for i = 1:size(ax,1)
            for j = 1:size(ax,2)
                iter = iter + 1;
                if iter > num, break, end
    
                cax = get_axis(fig,ax{i,j});
                cla(cax)
                hold on
                ydata = s.(fields_1D{l})(k).(q{iter});
                plot(s.(fields_1D{l})(k).r,ydata,'.')
                plot([s.(fields_hole{l})(k).p_xL s.(fields_hole{l})(k).p_xL],[min(ydata) max(ydata)])
                plot([s.(fields_hole{l})(k).p_xR s.(fields_hole{l})(k).p_xR],[min(ydata) max(ydata)])
                cax.PlotBoxAspectRatio = [1 1 1];
                if i == size(ax,1), xlabel('x (cm)'), end
                title(qstr{iter},'FontWeight','normal')
                an.String = ['t = ' num2str(s.(fields_1D{l})(k).t*1e6) '\mus = ' num2str(s.(fields_1D{l})(k).t/s.tau_hole,'%.2g') ' \tau'];
            end
        end
        frames{k} = getframe(fig);
    end
    close(fig)
    write_video([s.folder filesep fields_hole{l} '-transects'],frames);
end

%% Plot Hole Position vs Analytic Model

num = 1;
[fig,~,an] = open_subplot(num,'Visible','on');
an.Position = [0.1567    0.9611    0.7230    0.0339];
hold on
grid on
grid minor
if exp_file_found
    l = get_line_specs(6);
else
    l = get_line_specs(4);
end

lc = l(1);
plot([s.sim_hole.t].*1e6,[s.sim_hole.p_xL],'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
lc = l(2);
plot([s.sim_hole.t].*1e6,[s.sim_hole.p_xR],'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
xan = s.sim_int.hole([s.sim_hole.t]);
xan_l = xan + mean([s.sim_hole(1).p_xL s.sim_hole(1).p_xR]);
xan_r = -xan + mean([s.sim_hole(1).p_xL s.sim_hole(1).p_xR]);
lc = l(3);
plot([s.sim_hole.t].*1e6,xan_l,'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
lc = l(4);
plot([s.sim_hole.t].*1e6,xan_r,'LineWidth',2,'MarkerSize',2,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
if exp_file_found
    lc = l(5);
    plot([s.exp_hole.t].*1e6,[s.exp_hole.p_xL],'o','LineWidth',2,'MarkerSize',4,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
    lc = l(6);
    plot([s.exp_hole.t].*1e6,[s.exp_hole.p_xR],'o','LineWidth',2,'MarkerSize',4,'Color',lc.col,'MarkerFaceColor',lc.col,'MarkerEdgeColor',lc.col)
    xan = s.exp_int.hole([s.exp_hole.t]);
    xan_l = xan + mean([s.exp_hole(1).p_xL s.exp_hole(1).p_xR]);
    xan_r = -xan + mean([s.exp_hole(1).p_xL s.exp_hole(1).p_xR]);
end

if exp_file_found
    lgd = legend({'Sim: x_L','Sim: x_R','Model: x_L','Model: x_R','Exp: x_L','Exp: x_R'});
else
    lgd = legend({'Fit: x_L','Fit: x_R','Model: x_L','Model: x_R'});
end


lgd.Position = [0.7255    0.2929    0.2291    0.4417];
xlabel('t (\mus)')
ylabel('Hole Pos. (cm)')
saveas(fig,[s.folder filesep 'hole-prop.png'])
close(fig)

% save([s.folder filesep 'ion-hole-data.mat'],'s');

%%

% opening a subplot
num = length(t_plot);
[fig,ax,an] = open_subplot(num);
iter = 0;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end
        
        cax = get_axis(fig,ax{i,j});

        hold off
        l = get_line_specs(length(fields_1D));
        xdata = [s.(fields_1D{1})(closest_time_point(ind_t_plot(iter))).r];
        ydata = [s.(fields_1D{1})(closest_time_point(ind_t_plot(iter))).n];
        ydata = ydata./max(ydata);
        plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
        hold on
        xdata = [s.(fields_1D{2})(ind_t_plot(iter)).r];
        ydata = [s.(fields_1D{2})(ind_t_plot(iter)).n];
        ydata = ydata./max(ydata);
        plot(xdata,ydata,'LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
        hold on

        cax.FontSize = 11;
        if i == size(ax,1), xlabel('r (cm)'), end
        if j == 1, ylabel('n / n_0'), end
        title(['t = ' num2str(s.sim_1D(closest_time_point(ind_t_plot(iter))).t*1e6,'%.2g') ' \mus'],'FontWeight','normal')
        ylim([0 1.1].*max(ydata))
    end
end
an.String = [];

lgd = legend({'Sim','Exp'});
lgd.Position = [0.793153505527392,0.380132906369903,0.108849558323886,0.067314489047856];

end