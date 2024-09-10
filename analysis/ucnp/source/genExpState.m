function [] = genExpState(path,flags)
%% Set Up Directories and Other Program Options
set_path = [path filesep 'set_' num2str(flags.set)]; % path to where simulation files will be saved
mkdir(set_path);
save([setpath filesep 'exp-flags.mat'],'flags','-mat');
state_path = [set_path filesep 'init.txt']; % extension to be updated to .state later
config_path = [set_path filesep 'ucnp.txt']; % extension to be updated to .config later
settings_path = [set_path filesep 'plasma.txt']; % extension to be updated to .settings later

%% Load Data
config = readConfigFile(flags.config_path);
settings = readSettingsFile(flags.settings_path);
s = loadExpData(path,'os.mat');

ind = s.map(1).R < .01;
Ti_from_lif = mean(s.map(1).Ti(ind));
if ~isempty(flags.Ti), Ti_from_lif = flags.Ti; end

%% Trim Plasma Images

for i = 1:length([s.img.t])
    indx = s.img(i).x < flags.trim_exp_domain(1) & s.img(i).x > -flags.trim_exp_domain(1);
    indy = s.img(i).y < flags.trim_exp_domain(2) & s.img(i).y > -flags.trim_exp_domain(2);
    s.img(i).x = s.img(i).x(indx);
    s.img(i).y = s.img(i).y(indy);
    s.img(i).n = s.img(i).n(indy,indx);
    s.img(i).n_x = s.img(i).n_x(indx);
    s.img(i).n_y = s.img(i).n_y(indy);
    s.img(i).n_sg = s.img(i).n_sg(indy,indx);
    s.img(i).n_x_sg = s.img(i).n_x_sg(indx);
    s.img(i).n_y_sg = s.img(i).n_y_sg(indy);
end

%% Obtain State Background with Fit to Image
% fit density distribution with analytic distribution
density_min = max(s.map(1).n_hpr,[],'all')*flags.n_min;
if strcmp(flags.n_dist,'gauss')
    [n_fit] = fitImgWithGaussian(s.map(1).x_hpr,s.map(1).y_hpr,s.map(1).n_hpr,density_min);
elseif strcmp(flags.n_dist,'exp')
    [n_fit] = fitImgWithExp(s.map(1).x_hpr,s.map(1).y_hpr,s.map(1).n_hpr,density_min);
end

% plot fit to density distribution
fields = {'img','imgfit','imgres'};
fields_text = {'Img','Fit','Res'};

num = length(fields);
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];

ind = strcmp(fields,'imgfit');
clim = max(n_fit.(fields{ind}),[],'all');

iter = 0;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end
        cax = get_axis(fig,ax{i,j});
        xdata = n_fit.x;
        ydata = n_fit.y;
        zdata = n_fit.(fields{iter});
        if j == 2, clim = max(zdata,[],'all'); end
        imagesc(xdata,ydata,zdata)
        colorbar
        cax.YDir = 'Normal';
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 10;
        cax.CLim = [0 clim];
        if strcmp(fields{iter},'imgres')
            cax.CLim = [-clim clim]./10;
        end
        if i == size(ax,1), xlabel('x (cm)'), end
        if j == 1, ylabel('y (cm)'), end
        title(fields_text{iter},'FontWeight','normal')
    end
end
an.String = 'Fit to Initial Density Distribution';
exportgraphics(fig,[set_path filesep 'bgd-fit.png'],'Resolution',300);
close(fig)

%% Extract Te From Size Evolution

% extract sizes from fit
fit = struct;
fit(length([s.img.t])).sig_x = [];
fit(length([s.img.t])).sig_y = [];
for i = 1:length([s.img.t])
    if strcmp(flags.n_dist,'gauss')
    [results] = fitImgWithGaussian(s.map(i).x_hpr,s.map(i).y_hpr,s.map(i).n_hpr,density_min);
elseif strcmp(flags.n_dist,'exp')
    [results] = fitImgWithExp(s.map(i).x_hpr,s.map(i).y_hpr,s.map(i).n_hpr,density_min);
end
    fit(i).sig_x = results.sigx;
    fit(i).sig_y = results.sigy;
    fit(i).sig_2D = (fit(i).sig_x*fit(i).sig_y)^(1/2); % 1/3 because this is the experimental data
    fit(i).sig_3D = (fit(i).sig_x*fit(i).sig_y^2)^(1/3); % 1/3 because this is the experimental data
end

tau_2D = @(T) getTauExp(fit(1).sig_2D,T);
fitmodel = @(c,data) fit(1).sig_2D*sqrt(1+data.^2/tau_2D(c)^2);
p0 = settings.Te;
data = [s.img.t];
data = data';
zdata = [fit.sig_2D]';
T_fit = lsqcurvefit(fitmodel,p0,data,zdata,0,200);
Te_fit = T_fit - Ti_from_lif;

num = 1;
[fig,~,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1479    0.8056    0.5164    0.1020];
hold on
l = get_line_specs(2);  
plot(data./tau_2D(T_fit),fitmodel(T_fit,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data./tau_2D(T_fit),[fit.sig_2D],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
xlabel('t / \tau')
ylabel('\sigma (cm)')
an.String = ['T_e = ' num2str(T_fit) ' K'];
grid on 
grid minor
exportgraphics(fig,[path filesep 'sig2DvsTime.png'],'Resolution',300);
close(fig)

% decide which electron temperature to use
if isempty(flags.Te)
    if strcmp(flags.n_dist,'gauss')
        Te_for_files = Te_fit;
    elseif strcmp(flags.n_dist,'exp')
        Te_for_files = s.Te0;
    else
        error('<n_dist> must be either <gauss> or <exp>.')
    end
else
    Te_for_files = flags.Te;
end

%% Interpolate Image onto Larger, Non-Uniform Domain and Apply SG Filter to Background
% interpolate data onto larger uniform grid, apply sg filter, then interpolate onto non-uniform domain if necessary

% generate domain positions
if flags.is_uniform
    grid_growth = 1;
    grid_spread = 1;
else
    grid_growth = flags.grid_growth;
    grid_spread = flags.grid_spread;
end
[state.x,state.dx] = getNonUniformGrids(flags.Nx,flags.sim_domain(1),grid_growth,grid_spread);
[state.y,state.dy] = getNonUniformGrids(flags.Ny,flags.sim_domain(2),grid_growth,grid_spread);
[state.X,state.Y] = meshgrid(state.x,state.y);
[state.dX,state.dY] = meshgrid(state.dx,state.dy);
state.R = sqrt(state.X.^2+state.Y.^2);

state.n = interp2(s.map(1).x_hpr,s.map(1).y_hpr,s.map(1).n_hpr,state.X,state.Y) - 2.5e7;

n_bgd = n_fit.fit(state.X,state.Y);
for i = 1:numel(state.n)
    if state.R(i) > flags.bgd_radius || state.n(i) < density_min || isnan(state.n(i))
        state.n(i) = n_bgd(i);
    end
end

if flags.use_init_velocity
    state.mom_x = cts.cgs.mI*state.n*interp2(s.map(1).x_hpr,s.map(1).y_hpr,s.map(1).v_hpr,state.X,state.Y,'linear',0);
else
    state.mom_x = zeros(size(state.X));
end
for i = 1:numel(state.n)
    if state.R(i) > .35
        state.mom_x(i) = 0;
    end
end

%% Temperature Distribution

T_dist = state.n.^(1/1.5);
T_dist = T_dist./max(T_dist,[],'all');

%% Handle Settings File
% need to copy .settings path into set path, but first must determine the central value of n, Ti, and Te and the plasma size
% on each axis

% determine values to go into plasma.settings file
settings.n = n_fit.amp;
if strcmp(config.eq_set,'ideal_mhd') || strcmp(config.eq_set,'ideal_mhd_cons')
    settings.Ti = Te_for_files;
else
    settings.Ti = Ti_from_lif;
end
settings.Te = Te_for_files;
settings.sig_x = n_fit.sigx;
settings.sig_y = n_fit.sigy;
settings.x_lim = flags.sim_domain(1);
settings.y_lim = flags.sim_domain(2);
settings.Nx = flags.Nx;
settings.Ny = flags.Ny;
settings.n_min = density_min;

% write settings file
settings_data = cell(1e2,1);
fields = fieldnames(settings);
for i = 1:length(fields)
    settings_data{i} = [fields{i} ' = cgs = ' num2str(settings.(fields{i}))];
end
settings_data(i+1:end) = []; % remove excess cells
writecell(settings_data,settings_path,'QuoteStrings','none');
movefile(settings_path,[extractBefore(settings_path,'.txt') '.settings']);

%% Handle Config File
% need to update the duration and time_output_interval lines of the .config file to be consistent with the temporal of
% information of interest for the given data set

% the end time for the simulation is taken to be equal to the last experimental time point by default
total_temp = settings.Te + settings.Ti;
if isempty(flags.t_max)
    t_max = max(s.t)*sqrt(s.Te0/total_temp); % time in cgs
else
    t_max = flags.t_max*tau_2D(total_temp);
end
interval = t_max/flags.num_time_pts; % time interval between recording mhd grids

% read config file and replace the two aforementioned lines
C = readfile(flags.config_path);
ind = find(startsWith(C,'duration'));
if length(ind) ~= 1, error('duration line not found'); end
C{ind} = ['duration = ' num2str(t_max)];
ind = find(startsWith(C,'time_output_interval'));
if length(ind) ~= 1, error('time_output_interval line not found'); end
C{ind} = ['time_output_interval = ' num2str(interval)];

% write .config file
writecell(C,config_path,'QuoteStrings','none');
movefile(config_path,[extractBefore(config_path,'.txt') '.config']);

%% Write State File

% plot initial density profile
[fig,ax,an] = open_subplot(1,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
Nx_uni = (max(state.x)-min(state.x))/min(diff(state.x));
Ny_uni = (max(state.y)-min(state.y))/min(diff(state.y));
x_uni = linspace(min(state.x),max(state.x),Nx_uni);
y_uni = linspace(min(state.y),max(state.y),Ny_uni);
[X_uni,Y_uni] = meshgrid(x_uni,y_uni);
uni_grid = interp2(state.x,state.y,state.n,X_uni,Y_uni);
imagesc(x_uni,y_uni,uni_grid./1e8)
cax = ax{1};
cax.YDir = 'Normal';
cax.PlotBoxAspectRatio = [1 1 1];
cax.FontSize = 10;
colorbar
xlabel('x (cm)')
ylabel('y (cm)')
title('Initial Density (10^8 cm^-^3)','FontWeight','normal')
exportgraphics(fig,[set_path filesep 'init-n.png'],'Resolution',300);
close(fig)

% vars that need to be generated
grids.d_x = state.dX;
grids.d_y = state.dY;
grids.pos_x = state.X;
grids.pos_y = state.Y;
grids.be_x = zeros(size(state.X));
grids.be_y = zeros(size(state.X));
grids.be_z = zeros(size(state.X));
if strcmp(config.eq_set,'ideal_mhd')
    grids.rho = settings.m_i.*state.n;
    grids.temp = settings.Te*ones(size(state.X));
    grids.mom_x = state.mom_x;
    grids.mom_y = zeros(size(state.X));
    grids.mom_z = zeros(size(state.X));
    grids.bi_x = zeros(size(state.X));
    grids.bi_y = zeros(size(state.X));
    grids.bi_z = zeros(size(state.X));
    grids.grav_x = zeros(size(state.X));
    grids.grav_y = zeros(size(state.X));
elseif strcmp(config.eq_set,'ideal_2F')
    grids.i_rho = settings.m_i.*state.n;
    grids.e_rho = cts.cgs.mE.*state.n;
    grids.i_mom_x = state.mom_x;
    grids.i_mom_y = zeros(size(state.X));
    grids.e_mom_x = state.mom_x;
    grids.e_mom_y = zeros(size(state.X));
    grids.i_temp = settings.Ti*T_dist;
    grids.e_temp = settings.Te*ones(size(state.X));
    grids.bi_x = zeros(size(state.X));
    grids.bi_y = zeros(size(state.X));
    grids.bi_z = zeros(size(state.X));
    grids.E_x = zeros(size(state.X));
    grids.E_y = zeros(size(state.X));
    grids.E_z = zeros(size(state.X));
    grids.grav_x = zeros(size(state.X));
    grids.grav_y = zeros(size(state.X));
elseif strcmp(config.eq_set,'ideal_mhd_2E')
    grids.rho = settings.m_i.*state.n;
    grids.i_temp = settings.Ti*T_dist;
    grids.e_temp = settings.Te*ones(size(state.X));
    grids.mom_x = state.mom_x;
    grids.mom_y = zeros(size(state.X));
    grids.bi_x = zeros(size(state.X));
    grids.bi_y = zeros(size(state.X));
    grids.grav_x = zeros(size(state.X));
    grids.grav_y = zeros(size(state.X));
end

% begin writing state file
state_data = cell(1e5,1);
state_data{1} = 'xdim,ydim';
state_data{2} = [num2str(settings.Nx) ',' num2str(settings.Ny)];
state_data{3} = 'ion_mass';
state_data{4} = num2str(settings.m_i);
state_data{5} = 'adiabatic_index';
state_data{6} = num2str(settings.adiabatic_index);
state_data{7} = 't=0';
% write grids to state file
fields = fieldnames(grids);
iter = 7;
for i = 1:length(fields)
    iter = iter + 1;
    state_data{iter} = fields{i};
    for j = 1:size(grids.(fields{i}),2)
        iter = iter + 1;
        allOneString = sprintf('%.15g,', grids.(fields{i})(:,j)');
        allOneString = allOneString(1:end-1);% strip final comma
        state_data{iter} = allOneString;
    end
end
% remove excess cells inside state data
for i = 1:length(state_data)
    if isempty(state_data{i})
        state_data = state_data(1:i-1);
        break;
    end
end
% write file
writecell(state_data,state_path,'QuoteStrings','none');
movefile(state_path,[extractBefore(state_path,'.txt') '.state'],'f');


end