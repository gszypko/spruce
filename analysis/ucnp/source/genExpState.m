function [] = genExpState(path,flags)
%% Set Up Directories and Other Program Options
f = filesep;
set_path = [path f 'set_' num2str(flags.set)]; % path to where simulation files will be saved
mkdir(set_path);
state_path = [set_path f 'init.txt']; % extension to be updated to .state later
config_path = [set_path f 'ucnp.txt']; % extension to be updated to .config later
settings_path = [set_path f 'plasma.txt']; % extension to be updated to .settings later

%% Load Data
% load config and settings files from MHD GitHub, as specified by the paths within 'flags'
config = readConfigFile(flags.config_path);
settings = readSettingsFile(flags.settings_path);
eq_set = config.eq_set;
% load the experimental data and identify the smallest time point, which is to be used as the initial time point
load([path f 'os.mat'],'os');
[~,ind_t] = sort([os.delays]); 
os = os(ind_t);

%% Filter Plasma Images and Trim Domain
% retrieve images and domain information, trim domain, and apply 
imgs(length([os.delays])).t = [];
for i = 1:length([os.delays])
    % retrieve information from output structure
    imgs(i).t = os(i).delays*1e-9;
    imgs(i).x = os(i).imgs.xRelInMM./10;
    imgs(i).y = os(i).imgs.yRelInMM./10;
    imgs(i).n = os(i).imgs.density.*1e8;
    
    % determine number of kernels to be used with SG filter
    dx = mean(diff(imgs(i).x));
    dy = mean(diff(imgs(i).y));
    sg_kx = ceil(flags.sg_filt_length/dx/2)*2+1;
    sg_ky = ceil(flags.sg_filt_length/dy/2)*2+1;
    if sg_kx < 5, sg_kx = 5; end
    if sg_ky < 5, sg_ky = 5; end
    sg_bx = floor(sg_kx/2);
    sg_by = floor(sg_ky/2);

    % use sg filter on density images
    imgs(i).n_filt = sgfilt2D(imgs(i).n,sg_kx,sg_ky,3,3);
    imgs(i).n_filt_res = imgs(i).n_filt - imgs(i).n;

    % remove boundary cells of sg filter from os.imgs
    imgs(i).x = imgs(i).x(1+sg_bx:end-sg_bx);
    imgs(i).y = imgs(i).y(1+sg_by:end-sg_by);
    imgs(i).n = imgs(i).n(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    imgs(i).n_filt = imgs(i).n_filt(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    imgs(i).n_filt_res = imgs(i).n_filt_res(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);

    % trim domain
    indx = imgs(i).x < flags.trim_exp_domain(1) & imgs(i).x > -flags.trim_exp_domain(1);
    indy = imgs(i).y < flags.trim_exp_domain(2) & imgs(i).y > -flags.trim_exp_domain(2);
    imgs(i).x = imgs(i).x(indx);
    imgs(i).y = imgs(i).y(indy);
    imgs(i).n = imgs(i).n(indy,indx);
    imgs(i).n_filt = imgs(i).n_filt(indy,indx);
    imgs(i).n_filt_res = imgs(i).n_filt_res(indy,indx);
end

%% Plot Plasma Images
% open figure with the given panel number
num = 3;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
fields = {'n','n_filt','n_filt_res'};
for k = 1:length([imgs.t])
    disp(['Plotting: ' num2str(k) '/' num2str(length([imgs.t]))])
    % update axis information
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end

            cax = get_axis(fig,ax{i,j});
            xdata = imgs(k).x;
            ydata = imgs(k).y;
            zdata = imgs(k).(fields{j})./1e8;
            imagesc(xdata,ydata,zdata)
            colorbar
            cax.YDir = 'Normal';
            cax.PlotBoxAspectRatio = [1 1 1];
            cax.FontSize = 10;
            if i == size(ax,1), xlabel('x (cm)'), end
            if j == 1, ylabel('y (cm)'), end
            title('Density (10^8 cm^-^3)','FontWeight','normal')
            str1 = ['t = ' num2str(os(k).delays/1000,'%.2g') '\mus'];
            an.String = str1;
            an.Position = [0.196100917431193,0.872420262664165,0.589521675636846,0.116079737335835];
        end
    end
    % ensure limits for all subplots are set by filtered density distribution
    ind = strcmp(fields,'n');
    ind2 = strcmp(fields,'n_filt');
    ind3 = strcmp(fields,'n_filt_res');
    ax{ind}.CLim = ax{ind2}.CLim;
    ax{ind3}.CLim = [-1 1].*ax{ind2}.CLim(2)/2;
    % save figure
    saveas(fig,[path f 'imgs-' num2str(k) '.png']);
end
close(fig)

%% Obtain State Background with Fit to Image
if strcmp(flags.n_dist,'gauss')
    [n_fit] = fitImgWithGaussian(imgs(1).x,imgs(1).y,imgs(1).n,flags.n_min);
end

% plot fit to density distribution
fields = {'img','imgfit','imgres'};
fields_text = {'Raw','Fit','Res'};
num = length(fields);
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];

% update axis information
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
        if j == 3, cax.CLim = [-clim/10 clim/10]; end
        if i == size(ax,1), xlabel('x (cm)'), end
        if j == 1, ylabel('y (cm)'), end
        title(fields_text{iter},'FontWeight','normal')
    end
end
close(fig)

%% Extract Te From Size Evolution

% get central ion temperature from LIF fits
[~,ind_t] = min([os.delays]);
x = os(ind_t).map.x;
y = os(ind_t).map.y;
[X,Y] = meshgrid(x,y);
r = sqrt((X - n_fit.x0).^2+(Y - n_fit.y0).^2);
sig = (n_fit.sigx*n_fit.sigy)^(1/2);
ind_c = r < sig/5;
Ti_from_lif = mean(os(ind_t).map.Ti(ind_c));

% extract sizes from fit
fit = struct;
fit(length([imgs.t])).sig_x = [];
fit(length([imgs.t])).sig_y = [];
for i = 1:length([imgs.t])
    results = fitImgWithGaussian(imgs(i).x,imgs(i).y,imgs(i).n);
    fit(i).sig_x = results.sigx;
    fit(i).sig_y = results.sigy;
    fit(i).sig_2D = (fit(i).sig_x*fit(i).sig_y)^(1/2); % 1/3 because this is the experimental data
    fit(i).sig_3D = (fit(i).sig_x*fit(i).sig_y^2)^(1/3); % 1/3 because this is the experimental data
end

tau_2D = @(T) getTauExp(fit(1).sig_2D,T);
tau_3D = @(T) getTauExp(fit(1).sig_3D,T);
fitmodel = @(c,data) fit(1).sig_3D*sqrt(1+data.^2/tau_3D(c)^2);
p0 = settings.Te;
data = [imgs.t];
data = data';
zdata = [fit.sig_3D]';
T_fit = lsqcurvefit(fitmodel,p0,data,zdata,0,200);
Te_fit = T_fit - Ti_from_lif;

num = 1;
[fig,~,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1479    0.8056    0.5164    0.1020];
hold on
l = get_line_specs(2);  
plot(data./tau_3D(T_fit),fitmodel(T_fit,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data./tau_3D(T_fit),[fit.sig_3D],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
xlabel('t / \tau')
ylabel('\sigma (cm)')
an.String = ['T_e = ' num2str(Te_fit) ' K'];
grid on 
grid minor
saveas(fig,[path f 'Te.png']);
close(fig)

% decide which electron temperature to use
if isempty(flags.Te)
    Te_for_files = Te_fit;
else
    Te_for_files = flags.Te;
end

%% Interpolate Image onto Larger, Non-Uniform Domain and Apply SG Filter to Background
% interpolate data onto larger uniform grid, apply sg filter, then interpolate onto non-uniform domain if necessary

% create uniform domain for sg smoothing
x = linspace(-flags.sim_domain(1),flags.sim_domain(1),flags.num_grids);
y = linspace(-flags.sim_domain(2),flags.sim_domain(2),flags.num_grids);
dx = mean(diff(x)).*ones(size(x));
dy = mean(diff(y)).*ones(size(y));
[X,Y] = meshgrid(x,y);

% interpolate density onto grid
if flags.apply_sg_filt
    n_for_interpolation = imgs(1).n_filt;
else
    n_for_interpolation = imgs(1).n;
end
n = interp2(imgs(1).x,imgs(1).y,n_for_interpolation,X,Y);
n_bgd = n_fit.fit(X,Y);
for i = 1:length(y)
    for j = 1:length(x)
        if isnan(n(i,j)) || n(i,j) < flags.n_min
            n(i,j) = n_bgd(i,j);
        end
    end
end

% determine number of kernels to be used with SG filter
dx = mean(dx);
dy = mean(dy);
sg_kx = ceil(flags.sg_bgd_length/dx/2)*2+1;
sg_ky = ceil(flags.sg_bgd_length/dy/2)*2+1;
if sg_kx < 5, sg_kx = 5; end
if sg_ky < 5, sg_ky = 5; end
sg_bx = floor(sg_kx/2);
sg_by = floor(sg_ky/2);

% use sg filter on density images
n_sg = sgfilt2D(n,sg_kx,sg_ky,3,3);
for i = 1:size(n_sg,1)
    for j = 1:size(n_sg,2)
        if n_sg(i,j) < flags.n_min
            n_sg(i,j) = flags.n_min;
        end
    end
end

% generate domain positions
if flags.is_uniform
    state.x = x;
    state.y = y;
    state.dx = dx;
    state.dy = dy;
else
    [state.x,state.dx] = getNonUniformGrids(flags.num_grids,flags.sim_domain(1),flags.grid_growth,flags.grid_spread);
    [state.y,state.dy] = getNonUniformGrids(flags.num_grids,flags.sim_domain(2),flags.grid_growth,flags.grid_spread);
end
[state.X,state.Y] = meshgrid(state.x,state.y);
[state.dX,state.dY] = meshgrid(state.dx,state.dy);
state.R = sqrt(state.X.^2+state.Y.^2);
state.n = interp2(x,y,n,state.X,state.Y);

if flags.apply_sg_bgd
    for i = 1:size(n_sg,1)
        for j = 1:size(n_sg,2)
            if state.R(i,j) > flags.sg_radius
                state.n(i,j) = n_bgd(i,j);
            end
        end
    end
end

%% Handle Settings File
% need to copy .settings path into set path, but first must determine the central value of n, Ti, and Te and the plasma size
% on each axis

% determine values to go into plasma.settings file
settings.n = n_fit.amp;
if strcmp(eq_set,'ideal_mhd') || strcmp(eq_set,'ideal_mhd_cons')
    settings.Ti = Te_for_files;
else
    settings.Ti = Ti_from_lif;
end
settings.Te = Te_for_files;
settings.sig_x = n_fit.sigx;
settings.sig_y = n_fit.sigy;
settings.x_lim = flags.sim_domain(1);
settings.y_lim = flags.sim_domain(2);
settings.Nx = flags.num_grids;
settings.Ny = flags.num_grids;
settings.n_min = flags.n_min;

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
    t_max = max([imgs.t])*sqrt(os(1).Te/total_temp); % time in cgs
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
xlabel('y (cm)')
title('Density (10^8 cm^-^3)','FontWeight','normal')
saveas(fig,[set_path f 'init-n.png']);
close(fig)

% vars that need to be generated
grids.d_x = state.dX;
grids.d_y = state.dY;
grids.pos_x = state.X;
grids.pos_y = state.Y;
grids.be_x = zeros(size(state.X));
grids.be_y = zeros(size(state.X));
grids.be_z = zeros(size(state.X));
if strcmp(eq_set,'ideal_mhd')
    grids.rho = settings.m_i.*state.n;
    grids.temp = settings.Te*ones(size(state.X));
    grids.mom_x = zeros(size(state.X));
    grids.mom_y = zeros(size(state.X));
    grids.mom_z = zeros(size(state.X));
    grids.bi_x = zeros(size(state.X));
    grids.bi_y = zeros(size(state.X));
    grids.bi_z = zeros(size(state.X));
    grids.grav_x = zeros(size(state.X));
    grids.grav_y = zeros(size(state.X));
elseif strcmp(eq_set,'ideal_2F')
    grids.i_rho = settings.m_i.*state.n;
    grids.e_rho = cts.cgs.mE.*state.n;
    grids.i_mom_x = zeros(size(state.X));
    grids.i_mom_y = zeros(size(state.X));
    grids.e_mom_x = zeros(size(state.X));
    grids.e_mom_y = zeros(size(state.X));
    grids.i_temp = settings.Ti*ones(size(state.X));
    grids.e_temp = settings.Te*ones(size(state.X));
    grids.bi_x = zeros(size(state.X));
    grids.bi_y = zeros(size(state.X));
    grids.bi_z = zeros(size(state.X));
    grids.E_x = zeros(size(state.X));
    grids.E_y = zeros(size(state.X));
    grids.E_z = zeros(size(state.X));
    grids.grav_x = zeros(size(state.X));
    grids.grav_y = zeros(size(state.X));
elseif strcmp(eq_set,'ideal_mhd_2E')
    grids.rho = settings.m_i.*state.n;
    grids.i_temp = settings.Ti*ones(size(state.X));
    grids.e_temp = settings.Te*ones(size(state.X));
    grids.mom_x = zeros(size(state.X));
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