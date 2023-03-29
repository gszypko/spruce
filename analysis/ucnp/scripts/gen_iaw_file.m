%% User Controls
flags.Te = 48;
flags.sig_y = 0;
flags.set = 0;
flags.Nx = 401;
flags.Ny = 11;
flags.t_max = 39.3e-6;
flags.num_time_pts = 100;
flags.config_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.config';
flags.settings_path = 'C:\Users\grant\Documents\GitHub\mhd\ucnp.settings';
flags.n_iaw_sig = .2;

%% Set Up Directories and Other Program Options
main_path = 'C:\Users\grant\OneDrive\Research\mhd\projects\iaw-density-dist';
file_name = 'Killian 2012 - Fig 4d.csv';
file_name_no_space = file_name(~isspace(file_name));
folder_path = [main_path filesep extractBefore(file_name_no_space,'.csv')];
mkdir(folder_path)
set_path = [folder_path filesep 'set_' num2str(flags.set)]; % path to where simulation files will be saved
mkdir(set_path);
state_path = [set_path filesep 'init.txt']; % extension to be updated to .state later
config_path = [set_path filesep 'ucnp.txt']; % extension to be updated to .config later
settings_path = [set_path filesep 'plasma.txt']; % extension to be updated to .settings later
git_path = 'C:\Users\grant\Documents\GitHub\mhd';

% load config and settings files from MHD GitHub, as specified by the paths within 'flags'
config = readConfigFile(flags.config_path);
settings = readSettingsFile(flags.settings_path);
eq_set = config.eq_set;

%% Load Data

data = readcell([main_path filesep file_name]);
x_pos = cell2mat(data(:,1))'./10; % converts mm to cm
density = cell2mat(data(:,2))'.*1e15./1e6; % accounts for 10^15 units on plot axis, converts from m^-3 to cm^-3

[x_pos,ind_x] = sort(x_pos);
density = density(ind_x);

for i = 1:length(x_pos)-1
    if x_pos(i) == x_pos(i+1)
        density(i) = (density(i) + density(i+1))/2;
        x_pos(i+1) = [];
        density(i+1) = [];
    end
    if i > length(x_pos)-2, break; end
end

density_filt = sgolayfilt(density,3,9);

%% Fit Gaussian to Data

[p,n_filt_fit,n_filt_guess] = fitGaussian1D(x_pos,density_filt);

fig = figure;
plot(x_pos,density_filt)
hold on
plot(x_pos,n_filt_guess)
plot(x_pos,n_filt_fit)
legend({'data','guess','fit'})
close(fig)

%% Handle Settings File
% need to copy .settings path into set path, but first must determine the central value of n, Ti, and Te and the plasma size
% on each axis

% determine values to go into plasma.settings file
settings.n = p(1);
if strcmp(eq_set,'ideal_mhd') || strcmp(eq_set,'ideal_mhd_cons')
    settings.Ti = flags.Te;
else
    settings.Ti = 1;
end
settings.Te = flags.Te;
settings.sig_x = p(3);
settings.sig_y = 1e10;
settings.x_lim = max(x_pos);
settings.y_lim = 1;
settings.Nx = flags.Nx;
settings.Ny = flags.Ny;
settings.n_min = min(density);
settings.n_iaw_sig = flags.n_iaw_sig;

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
t_max = flags.t_max;
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

[state.x,state.dx] = getNonUniformGrids(flags.Nx,max(abs(x_pos)),1,1);
[state.y,state.dy] = getNonUniformGrids(flags.Ny,1,1,1);
[state.X,state.Y] = meshgrid(state.x,state.y);
[state.dX,state.dY] = meshgrid(state.dx,state.dy);
state.n = zeros(size(state.X));
for i = 1:size(state.n,1)
    state.n(i,:) = interp1(x_pos,density_filt,state.x,'makima');
end

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