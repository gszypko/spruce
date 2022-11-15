function [] = genExpState(path,flags)
%% Set Up Directories and Other Program Options
% generate set paths
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
[~,ind_t] = min([os.delays]); 

%% Filter Out Data Far From Plasma Center
for i = 1:length([os.delays])
    indx = os(i).imgs.xRelInMM < flags.exp_window(1)*10 & os(i).imgs.xRelInMM > -flags.exp_window(1)*10;
    indy = os(i).imgs.yRelInMM < flags.exp_window(2)*10 & os(i).imgs.yRelInMM > -flags.exp_window(2)*10;
    os(i).imgs.xRelInMM = os(i).imgs.xRelInMM(indx);
    os(i).imgs.yRelInMM = os(i).imgs.yRelInMM(indy);
    os(i).imgs.density = os(i).imgs.density(indy,indx);
end

%% Apply SG Filter to Integrated Density and Filter Hot Pixels
% for each time point, 
for i = 1:length([os.delays])
    % determine number of kernels to be used with SG filter
    dx = mean(diff(os(i).imgs.xRelInMM))/10;
    dy = mean(diff(os(i).imgs.yRelInMM))/10;
    sg_kx = ceil(flags.sg_imgs_length/dx/2)*2+1;
    sg_ky = ceil(flags.sg_imgs_length/dy/2)*2+1;
    if sg_kx < 5, sg_kx = 5; end
    if sg_ky < 5, sg_ky = 5; end
    sg_bx = floor(sg_kx/2);
    sg_by = floor(sg_ky/2);
    % use sg filter on density images
    os(i).imgs.density_filt = sgfilt2D(os(i).imgs.density,sg_kx,sg_ky,3,3);
    os(i).imgs.density_filt_res = os(i).imgs.density_filt - os(i).imgs.density;
    % remove boundary cells of sg filter from os.imgs
    os(i).imgs.xRelInMM = os(i).imgs.xRelInMM(1+sg_bx:end-sg_bx);
    os(i).imgs.yRelInMM = os(i).imgs.yRelInMM(1+sg_by:end-sg_by);
    os(i).imgs.density = os(i).imgs.density(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    os(i).imgs.density_filt = os(i).imgs.density_filt(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    os(i).imgs.density_filt_res = os(i).imgs.density_filt_res(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    % bin image if it has too many pixels
    fields = {'density','density_filt','density_filt_res'};
    for j = 1:length(fields)
        [xbin,ybin,os(i).imgs.(fields{j})] = binImgForFit(os(i).imgs.xRelInMM,os(i).imgs.yRelInMM,os(i).imgs.(fields{j}),600);
    end
    os(i).imgs.xRelInMM = xbin;
    os(i).imgs.yRelInMM = ybin;
end

%% Plot Density Image Derived from Spectrum Integration
% open figure with the given panel number
num = 3;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
fields = {'density','density_filt','density_filt_res'};
for k = 1:length(os)
    disp(['Plotting: ' num2str(k) '/' num2str(length(os))])
    % update axis information
    iter = 0;
    for i = 1:size(ax,1)
        for j = 1:size(ax,2)
            iter = iter + 1;
            if iter > num, break, end

            cax = get_axis(fig,ax{i,j});
            xdata = os(k).imgs.xRelInMM./10;
            ydata = os(k).imgs.yRelInMM./10;
            zdata = os(k).imgs.(fields{j});
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
    ind = strcmp(fields,'density');
    ind2 = strcmp(fields,'density_filt');
    ind3 = strcmp(fields,'density_filt_res');
    ax{ind}.CLim = ax{ind2}.CLim;
    ax{ind3}.CLim = [-1 1].*ax{ind2}.CLim(2)/2;
    % save figure
    saveas(fig,[path f 'imgs-' num2str(k) '.png']);
end
close(fig)

%% Get Density from LIF Fits
% for each time point, 
for i = 1:length([os.delays])
    % determine number of kernels to be used with SG filter
    dx = mean(diff(os(i).map.x))/10;
    dy = mean(diff(os(i).map.y))/10;
    sg_kx = ceil(flags.sg_imgs_length/dx/2)*2+1;
    sg_ky = ceil(flags.sg_imgs_length/dy/2)*2+1;
    if sg_kx < 5, sg_kx = 5; end
    if sg_ky < 5, sg_ky = 5; end
    sg_bx = floor(sg_kx/2);
    sg_by = floor(sg_ky/2);
    % use sg filter on density images
    os(i).map.nFit_filt = sgfilt2D(os(i).map.nFit,sg_kx,sg_ky,3,3);
    os(i).map.nFit_res = os(i).map.nFit_filt - os(i).map.nFit;
    % remove boundary cells of sg filter from os.map
    os(i).map.x = os(i).map.x(1+sg_bx:end-sg_bx);
    os(i).map.y = os(i).map.y(1+sg_by:end-sg_by);
    os(i).map.nFit = os(i).map.nFit(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    os(i).map.nFit_filt = os(i).map.nFit_filt(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    os(i).map.nFit_res = os(i).map.nFit_res(1+sg_bx:end-sg_bx,1+sg_by:end-sg_by);
    % bin image if it has too many pixels
    fields = {'nFit','nFit_filt','nFit_res'};
    for j = 1:length(fields)
        [xbin,ybin,os(i).map.(fields{j})] = binImgForFit(os(i).map.x,os(i).map.y,os(i).map.(fields{j}),600);
    end
    os(i).map.x = xbin;
    os(i).map.y = ybin;
end

%% Extract Te From Size Evolution

% for each time point, 
fit = struct;
fit(length([os.delays])).sig_x = [];
fit(length([os.delays])).sig_y = [];
for k = 1:length([os.delays])
    x = os(k).imgs.xRelInMM/10;
    y = os(k).imgs.yRelInMM/10;
    img = os(k).imgs.density;
    results = fitImgWithGaussian(x,y,img);
    fit(k).sig_x = results.sigx;
    fit(k).sig_y = results.sigy;
    fit(k).sig = (fit(k).sig_x*fit(k).sig_y^2)^(1/3);
    fit(k).sig_norm = fit(1).sig;
end

tau = @(Te) sqrt(cts.cgs.mI*fit(1).sig^2/cts.cgs.kB/Te);
fitmodel = @(c,data) fit(1).sig*sqrt(1+data.^2/tau(c)^2);
p0 = settings.Te;
data = [os.delays].*1e-9;
data = data';
zdata = [fit.sig]';
[p,~,R,~,~,~,J] = lsqcurvefit(fitmodel,p0,data,zdata);

num = 1;
[fig,ax,an] = open_subplot(num,'Visible',flags.figvis);
an.Position = [0.1479    0.8056    0.5164    0.1020];
cax = ax{1};
hold on
l = get_line_specs(2);  
plot(data./tau(p),fitmodel(p,data),'LineWidth',2,'MarkerSize',4,'Color',l(1).col,'MarkerFaceColor',l(1).col,'MarkerEdgeColor',l(1).col)
plot(data./tau(p),[fit.sig],'o','LineWidth',2,'MarkerSize',4,'Color',l(2).col,'MarkerFaceColor',l(2).col,'MarkerEdgeColor',l(2).col)
an.String = ['T_e = ' num2str(p) ' K'];
saveas(fig,[path f 'Te.png']);
close(fig)

%% Handle Config File
% need to update the duration and time_output_interval lines of the .config file to be consistent with the temporal of
% information of interest for the given data set

% the end time for the simulation is taken to be equal to the last experimental time point by default
t_max = max([os.delays])*1e-9*sqrt(os(1).Te/settings.Te); % time in cgs
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

%% Handle Settings File
% need to copy .settings path into set path, but first must determine the central value of n, Ti, and Te and the plasma size
% on each axis

% fit plasma with Gaussian to obtain its size
[fit] = fitImgWithGaussian(os(ind_t).imgs.xRelInMM/10,os(ind_t).imgs.yRelInMM/10,os(ind_t).imgs.density);

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
        xdata = fit.x;
        ydata = fit.y;
        zdata = fit.(fields{iter});
        imagesc(xdata,ydata,zdata)
        colorbar
        cax.YDir = 'Normal';
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 10;
        if i == size(ax,1), xlabel('x (cm)'), end
        if j == 1, ylabel('y (cm)'), end
        title(fields_text{iter},'FontWeight','normal')
    end
end
close(fig)

% find cells nearest the plasma center
x = os(ind_t).map.x;
y = os(ind_t).map.y;
[X,Y] = meshgrid(x,y);
r = sqrt((X - fit.x0).^2+(Y - fit.y0).^2);
sig = (fit.sigx*fit.sigy^2)^(1/3);
ind_c = r < sig/5;

% determine values to go into plasma.settings file
settings.n = fit.amp;
settings.Ti = mean(os(ind_t).map.Ti(ind_c));
if isempty(flags.Te)
    settings.Te = os.Te;
else
    settings.Te = flags.Te;
end
settings.sig_x = fit.sigx;
settings.sig_y = fit.sigy;
settings.x_lim = flags.sim_window(1);
settings.y_lim = flags.sim_window(2);
settings.Nx = flags.Nx;
settings.Ny = flags.Ny;
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

%% Write State File

% option for whether to use smoothed density distribution or not
if flags.smooth_density
    if ~flags.useLIFFits
        density = os(ind_t).imgs.density_filt;
    else
        density = os(ind_t).map.nFit_filt;
    end
else
    if ~flags.useLIFFits
        density = os(ind_t).imgs.density;
    else
        density = os(ind_t).map.nFit;
    end
end

% enforce density minimum on the density grid prior to doing any smoothing
for i = 1:size(density,1)
    for j = 1:size(density,2)
        if density(i,j) < settings.n_min/1e8
            density(i,j) = settings.n_min/1e8;
        end
    end
end

% identify spatial extent of experimental domain in cgs
x_L = min(os(ind_t).imgs.xRelInMM)/10; 
x_R = max(os(ind_t).imgs.xRelInMM)/10;
y_L = min(os(ind_t).imgs.yRelInMM)/10;
y_R = max(os(ind_t).imgs.yRelInMM)/10;

% determine if exp domain is smaller than sim domain
is_exp_domain_smaller = -settings.x_lim < x_L && settings.x_lim > x_R && -settings.y_lim < y_L && settings.y_lim > y_R;

% determine number of kernels to be used with SG filter
x = linspace(-settings.x_lim,settings.x_lim,settings.Nx);
y = linspace(-settings.y_lim,settings.y_lim,settings.Ny);
sg_length = flags.sg_back_length;
dx = mean(diff(x));
dy = mean(diff(y));
sg_kx = ceil(sg_length/dx/2)*2+1;
sg_ky = ceil(sg_length/dy/2)*2+1;
sg_bx = floor(sg_kx/2);
sg_by = floor(sg_ky/2);

% get spatial grids to account for sg_filter cutoff
state.x = linspace(-settings.x_lim-sg_bx*dx,settings.x_lim+sg_bx*dx,settings.Nx+2*sg_bx);
state.y = linspace(-settings.y_lim-sg_by*dy,settings.y_lim+sg_by*dy,settings.Ny+2*sg_by);
[X,Y] = meshgrid(state.x,state.y);
state.X = X;
state.Y = Y;

% interpolate experimental density distribution onto the simulation grid
if ~flags.useLIFFits
    X = os(ind_t).imgs.xRelInMM/10;
    Y = os(ind_t).imgs.yRelInMM/10;
else
    X = os(ind_t).map.x/10;
    Y = os(ind_t).map.y/10;
end
V = density*1e8;
Xq = state.X;
Yq = state.Y;
state.n = interp2(X,Y,V,Xq,Yq,'linear',settings.n_min);

% enforce density minimum again
for i = 1:size(state.n,1)
    for j = 1:size(state.n,2)
        if state.n(i,j) < settings.n_min
            state.n(i,j) = settings.n_min;
        end
    end
end

% apply sg filter to state grid and only keep the smoothed values for positions exterior to simulation domain or less than min
x_L_sg = x_L + sg_length/2; 
x_R_sg = x_R - sg_length/2;
y_L_sg = y_L + sg_length/2;
y_R_sg = y_R - sg_length/2;
state.n_sg = sgfilt2D(state.n,sg_kx,sg_ky,1,1);
for i = 1:size(state.n,1)
    for j = 1:size(state.n,2)
        is_in_sg_domain = state.X(i,j) < x_R_sg && state.X(i,j) > x_L_sg && state.Y(i,j) < y_R_sg && state.Y(i,j) > y_L_sg;
        is_above_min = state.n(i,j) > settings.n_min;
        if is_exp_domain_smaller
            if is_in_sg_domain && is_above_min
                state.n_sg(i,j) = state.n(i,j);
            end
        else
            if is_above_min
                state.n_sg(i,j) = state.n(i,j);
            end
        end
    end
end

% remove sg boundary cells
ind_x = 1+sg_bx:length(state.x)-sg_bx;
ind_y = 1+sg_by:length(state.y)-sg_by;
state.x = state.x(ind_x);
state.y = state.y(ind_y);
[X,Y] = meshgrid(state.x,state.y);
state.X = X;
state.Y = Y;
state.n = state.n(ind_y,ind_x);
state.n_sg = state.n_sg(ind_y,ind_x);

[fig,ax,an] = open_subplot(1,'Visible',flags.figvis);
an.Position = [0.1567    0.8909    0.7230    0.0801];
imagesc(state.x,state.y,state.n_sg)
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

% begin handling grids - variable names depends on equation set

% vars that need to be generated
grids.d_x = mean(diff(x))*ones(size(X));
grids.d_y = mean(diff(y))*ones(size(Y));
grids.pos_x = X;
grids.pos_y = Y;
grids.be_x = zeros(size(X));
grids.be_y = zeros(size(X));
if strcmp(eq_set,'ideal_mhd')
    grids.rho = settings.m_i.*state.n_sg;
    grids.temp = settings.Te*ones(size(X));
    grids.mom_x = zeros(size(X));
    grids.mom_y = zeros(size(X));
    grids.bi_x = zeros(size(X));
    grids.bi_y = zeros(size(X));
    grids.grav_x = zeros(size(X));
    grids.grav_y = zeros(size(X));
elseif strcmp(eq_set,'ideal_2F')
    grids.i_rho = settings.m_i.*state.n_sg;
    grids.e_rho = cts.cgs.mE.*state.n_sg;
    grids.i_mom_x = zeros(size(X));
    grids.i_mom_y = zeros(size(X));
    grids.e_mom_x = zeros(size(X));
    grids.e_mom_y = zeros(size(X));
    grids.i_temp = settings.Ti*ones(size(X));
    grids.e_temp = settings.Te*ones(size(X));
    grids.bi_x = zeros(size(X));
    grids.bi_y = zeros(size(X));
    grids.bi_z = zeros(size(X));
    grids.E_x = zeros(size(X));
    grids.E_y = zeros(size(X));
    grids.E_z = zeros(size(X));
    grids.grav_x = zeros(size(X));
    grids.grav_y = zeros(size(X));
elseif strcmp(eq_set,'ideal_mhd_2E')
    grids.rho = settings.m_i.*state.n_sg;
    grids.i_temp = settings.Ti*ones(size(X));
    grids.e_temp = settings.Te*ones(size(X));
    grids.mom_x = zeros(size(X));
    grids.mom_y = zeros(size(X));
    grids.bi_x = zeros(size(X));
    grids.bi_y = zeros(size(X));
    grids.grav_x = zeros(size(X));
    grids.grav_y = zeros(size(X));
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
        allOneString = sprintf('%.9g,', grids.(fields{i})(:,j)');
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