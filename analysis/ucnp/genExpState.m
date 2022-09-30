%% Script Initialization

% clean up workspace
clc, clearvars, close all, f = filesep; setpath;

% specify paths to relevant folders
target_path = 'C:\Users\grant\OneDrive\Research\mhd\data-expsims\an-mhd-09.26.22\HIH-4';
config_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.config';
settings_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.settings';

% generate set paths
set_id = 1;
set_path = [target_path f 'set_' num2str(set_id)];
mID = mkdir(set_path);
state_path_txt = [set_path f 'init.txt'];
state_path_state = [set_path f 'init.state'];
settings_path_txt = [set_path f 'plasma.txt'];
settings_path_settings = [set_path f 'plasma.settings'];
config_path_txt = [set_path f 'ucnp.txt'];
config_path_config = [set_path f 'ucnp.config'];

% image filter options
sg_kernel = 15; % number of grids to use for sg smoothing of density images
sg_kernel_lif = 7;
sg_kernel_state = 15;
gauss_kernel = 3; % sigma value for imguassfilt

% other script options
num_pts = 250;
plot_density_images = true;
plot_lif_fits = true;
plasma_distribution = 'gaussian'; % 'gaussian' or 'exponential'
figure_visibility = 'on'; % 'on' or 'off'

%% Load Data

config = readConfigFile(config_path);
eq_set = config.eq_set;
settings = readSettingsFile(settings_path);
load([target_path f 'os.mat']);
[~,ind_t] = min([os.delays]); % identify smallest time point, to be used for initialization

%% Smooth Data
sg_boundary = floor(sg_kernel/2);
sg_boundary_lif = floor(sg_kernel_lif/2);
for i = 1:length(os)
    % use sg filter on density images
    os(i).imgs.density_filt = sgfilt2D(os(i).imgs.density,sg_kernel,sg_kernel,3,3);
    os(i).imgs.density_filt_res = os(i).imgs.density_filt - os(i).imgs.density;
    % remove boundary cells of sg filter from os.imgs
    os(i).imgs.xRelInMM = os(i).imgs.xRelInMM(1+sg_boundary:end-sg_boundary);
    os(i).imgs.yRelInMM = os(i).imgs.yRelInMM(1+sg_boundary:end-sg_boundary);
    os(i).imgs.density = os(i).imgs.density(1+sg_boundary:end-sg_boundary,1+sg_boundary:end-sg_boundary);
    os(i).imgs.density_filt = os(i).imgs.density_filt(1+sg_boundary:end-sg_boundary,1+sg_boundary:end-sg_boundary);
    os(i).imgs.density_filt_res = os(i).imgs.density_filt_res(1+sg_boundary:end-sg_boundary,1+sg_boundary:end-sg_boundary);
    % use sg filter on lif fits
    fields = fieldnames(os(i).map);
    for j = 3:length(fields)
        os(i).map.(fields{j}) = sgfilt2D(os(i).map.(fields{j}),sg_kernel_lif,sg_kernel_lif,3,3);
        os(i).map.(fields{j}) = os(i).map.(fields{j})(1+sg_boundary_lif:end-sg_boundary_lif,1+sg_boundary_lif:end-sg_boundary_lif);
    end
    % truncate lif maps based on sg filter boundary cells
    os(i).map.x = os(i).map.x(1+sg_boundary_lif:end-sg_boundary_lif);
    os(i).map.y = os(i).map.y(1+sg_boundary_lif:end-sg_boundary_lif);

end

%% Plot Density Image Derived from Spectrum Integration
if plot_density_images
    % open figure
    panel = [1 3];
    num = panel(1)*panel(2);
    [fig,ax,an] = open_subplot(panel(1),panel(2),'Visible',figure_visibility,num);
    fig.Position = [7.916666666666666e+02,441,1.162666666666667e+03,3.553333333333333e+02];
    an.Position = [0.1595    0.9084    0.7230    0.0801];
    fields = {'density','density_filt','density_filt_res'};
    
    % plot frames and save video
    frames = cell(1,length(os));
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
                if strcmp(fields{j},'density_filt'), clims = cax.CLim; end
            end
        end
        for i = 1:size(ax,1)
            for j = 1:size(ax,2)
                cax = get_axis(fig,ax{i,j});
                cax.CLim = clims;
                if j == 3, cax.CLim = cax.CLim./10; end
            end
        end
        dlm = ' - ';
        frames{k} = getframe(fig);
        pause(1)
    end
    close(fig)
    write_video([target_path f 'density-from-integration-images'],frames);
end

%% Plot LIF Fits for Density, Velocity, and Temperature
if plot_lif_fits
    % open figure
    fields = {'nFit','vExp','Ti'};
    fields_title = {'Density (10^8 cm^-^3)','Velocity (cm/s)', 'Temperature (K)'};
    cut = 2;
    panel = [1 3];
    num = panel(1)*panel(2);
    [fig,ax,an] = open_subplot(panel(1),panel(2),'Visible',figure_visibility,num);
    fig.Position = [7.916666666666666e+02,441,1.162666666666667e+03,3.553333333333333e+02];
    an.Position = [0.1595    0.9084    0.7230    0.0801];
    
    % plot frames and save video
    frames = cell(1,length(os));
    for k = 1:length(os)
        disp(['Plotting: ' num2str(k) '/' num2str(length(os))])
        
        % update axis information
        iter = 0;
        for i = 1:size(ax,1)
            for j = 1:size(ax,2)
                iter = iter + 1;
                if iter > num, break, end
    
                cax = get_axis(fig,ax{i,j});
                xdata = os(k).map.x./10;
                ydata = os(k).map.y./10;
                zdata = os(k).map.(fields{j});
                if strcmp(fields{j},'vExp'), zdata = zdata./10; end
                imagesc(xdata(1+cut:end-cut),ydata(1+cut:end-cut),zdata(1+cut:end-cut,1+cut:end-cut))
                colorbar
                cax.YDir = 'Normal';
                cax.PlotBoxAspectRatio = [1 1 1];
                cax.FontSize = 10;
                if i == size(ax,1), xlabel('x (cm)'), end
                if j == 1, ylabel('y (cm)'), end
                title(fields_title{j},'FontWeight','normal')
                str1 = ['t = ' num2str(os(k).delays/1000,'%.2g') '\mus'];
                an.String = str1;
                an.Position = [0.196100917431193,0.872420262664165,0.589521675636846,0.116079737335835];
                if strcmp(fields{j},'density_filt'), clims = cax.CLim; end
            end
        end
        dlm = ' - ';
        frames{k} = getframe(fig);
        pause(1)
    end
    close(fig)
    write_video([target_path f 'lif-fits'],frames);
end



%% Handle Config File
% need to update the duration and time_output_interval lines of the .config file to be consistent with the temporal of
% information of interest for the given data set

% determine end time for simulation
t_max = max([os.delays])*1e-9;
interval = t_max/num_pts;

% read in config text
C = readfile(config_path);
ind = find(startsWith(C,'duration'));
if length(ind) ~= 1, error('duration line not found'); end
C{ind} = ['duration = ' num2str(t_max)];
ind = find(startsWith(C,'time_output_interval'));
if length(ind) ~= 1, error('time_output_interval line not found'); end
C{ind} = ['time_output_interval = ' num2str(interval)];

% write .config file
writecell(C,config_path_txt,'QuoteStrings','none');
movefile(config_path_txt,config_path_config);


%% Handle Settings File


% need to copy .settings path into set path, but first must determine the central value of n, Ti, and Te and the plasma size
% on each axis

% fit plasma with Gaussian to obtain its size
if strcmp(plasma_distribution,'gaussian')
    [fit] = fitImgWithGaussian(os(ind_t).imgs.xRelInMM,os(ind_t).imgs.yRelInMM,os(ind_t).imgs.density);
end

% plot fit to density distribution
fields = {'imgbin','imgfit','imgres'};
fields_text = {'Raw','Fit','Res'};
panel = [1 length(fields)];
num = length(fields);
[fig,ax,an] = open_subplot(panel(1),panel(2),'Visible',figure_visibility,num);
fig.Position = [1.082000000000000e+02,2.242000000000000e+02,1.162400000000000e+03,3.552000000000001e+02];
an.Position = [0.1595    0.9084    0.7230    0.0801];

% update axis information
iter = 0;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end

        cax = get_axis(fig,ax{i,j});
        xdata = fit.xRelInMM./10;
        ydata = fit.yRelInMM./10;
        zdata = fit.(fields{iter});
        imagesc(xdata,ydata,zdata)
        colorbar
        cax.YDir = 'Normal';
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 10;
        if i == size(ax,1), xlabel('x (cm)'), end
        if j == 1, ylabel('y (cm)'), end
        title(fields_text{iter},'FontWeight','normal')
        if strcmp(fields{j},'density_filt'), clims = cax.CLim; end
    end
end
pause(1)
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
settings.Te = os.Te;
settings.sig_x = fit.sigx/10;
settings.sig_y = fit.sigy/10;

% write settings file
settings_data = cell(1e2,1);
fields = fieldnames(settings);
for i = 1:length(fields)
    settings_data{i} = [fields{i} ' = cgs = ' num2str(settings.(fields{i}))];
end
settings_data(i+1:end) = []; % remove excess cells
writecell(settings_data,settings_path_txt,'QuoteStrings','none');
movefile(settings_path_txt,settings_path_settings);

%% Write State File

% generate spatial grids for density
sg_boundary_state = floor(sg_kernel_state/2);
x = linspace(-settings.x_lim,settings.x_lim,settings.Nx);
dx = mean(diff(x));
x_extended = linspace(-settings.x_lim-dx*sg_boundary_state,settings.x_lim+dx*sg_boundary_state,settings.Nx+sg_boundary_state*2);
y = linspace(-settings.y_lim,settings.y_lim,settings.Ny);
dy = mean(diff(x));
y_extended = linspace(-settings.y_lim-dy*sg_boundary_state,settings.y_lim+dy*sg_boundary_state,settings.Ny+sg_boundary_state*2);
[X,Y] = meshgrid(x_extended,y_extended);
R = sqrt(X.^2+Y.^2);
ind_r = R > 0.25;
density = settings.n_min*ones(size(X));
for i = 1:size(X,1)
    for j = 1:size(X,2)
        x = os(ind_t).imgs.xRelInMM./10;
        y = os(ind_t).imgs.yRelInMM./10;
        V = os(ind_t).imgs.density_filt*1e8;
        Xq = X(i,j);
        Yq = Y(i,j);
        if Xq < max(x) && Xq > min(x) && Yq < max(y) && Yq > min(y)
            value = interp2(x,y,V,Xq,Yq);
            if value > density(i,j), density(i,j) = value; end
        end
    end
end
density_filt = sgfilt2D(density,sg_kernel_state,sg_kernel_state,3,3);
for i = 1:size(density,1)
    for j = 1:size(density,2)
        if ind_r(i,j), density(i,j) = density_filt(i,j); end
    end
end
x = x_extended(1+sg_boundary_state:end-sg_boundary_state);
y = y_extended(1+sg_boundary_state:end-sg_boundary_state);
[X,Y] = meshgrid(x,y);
density = density(1+sg_boundary_state:end-sg_boundary_state,1+sg_boundary_state:end-sg_boundary_state);

imagesc(x,y,density)
ax = gca;
ax.YDir = 'Normal';

% begin handling grids - variable names depends on equation set

% vars that need to be generated
grids.d_x = mean(diff(x))*ones(size(X));
grids.d_y = mean(diff(y))*ones(size(Y));
grids.pos_x = X;
grids.pos_y = Y;
grids.be_x = zeros(size(X));
grids.be_y = zeros(size(X));
if strcmp(eq_set,'ideal_mhd')
    grids.rho = settings.m_i.*density;
    grids.temp = settings.Te*ones(size(X));
    grids.mom_x = zeros(size(X));
    grids.mom_y = zeros(size(X));
    grids.bi_x = zeros(size(X));
    grids.bi_y = zeros(size(X));
    grids.grav_x = zeros(size(X));
    grids.grav_y = zeros(size(X));
elseif strcmp(eq_set,'ideal_2F')
    grids.i_rho = settings.m_i.*density;
    grids.e_rho = settings.m_i.*density;
    grids.i_mom_x = zeros(size(X));
    grids.i_mom_y = zeros(size(X));
    grids.e_mom_x = zeros(size(X));
    grids.e_mom_y = zeros(size(X));
    grids.i_temp = settings.Ti*ones(size(X));
    grids.e_temp = settings.Te*ones(size(X));
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
        allOneString = sprintf('%.5g,', grids.(fields{i})(:,j)');
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
writecell(state_data,state_path_txt,'QuoteStrings','none');
movefile(state_path_txt,state_path_state);