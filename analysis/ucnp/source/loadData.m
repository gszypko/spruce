function [data] = loadData(folder,flags)
% load settings, config, and state files
data.folder = folder;
data.settings = readSettingsFile([folder filesep 'plasma.settings']);
data.config = readConfigFile([folder filesep 'ucnp.config']);
data.state.init = readStateFile([folder filesep 'init.state'],flags.ghost_cells,data.config.eq_set);

% determine adiabatic indices for each species
data.adiabatic_index_e = data.settings.adiabatic_index;
data.adiabatic_index_i = data.settings.adiabatic_index;
if data.config.global_temp_e, data.adiabatic_index_e = 1; end
if data.config.global_temp_i, data.adiabatic_index_i = 1; end

% determine grids to be loaded from mhd.out
if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
    grids_to_load = {'n','v_x','v_y','temp'};
elseif max(strcmp({'ideal_2F'},data.config.eq_set))
    grids_to_load = {'i_n','i_v_x','i_v_y','i_temp','e_temp'};
elseif strcmp(data.config.eq_set,{'ideal_mhd_2E'})
    grids_to_load = {'n','v_x','v_y','i_temp','e_temp'};
end
if flags.plot_grids
    grids_to_load = [grids_to_load flags.vars]; 
end
grids_to_load = unique(grids_to_load);

% load grids from mhd.out
data.grids = readOutFile([folder filesep 'mhd.out'],flags.ghost_cells,grids_to_load,flags.time_window,flags.time_interval);

% determine central electron and ion temperatures from settings file and compute tau_exp
data.Te = data.settings.Te;
if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
    data.Ti = data.settings.Te;
elseif max(strcmp({'ideal_2F','ideal_mhd_2E'},data.config.eq_set))
    data.Ti = data.settings.Ti;
else
    error('Error: equation set is not valid.')
end

% determine tau_exp from geometric mean of rms plasma size
data.sigx = data.settings.sig_x;
data.sigy = data.settings.sig_y;
data.sig0 = (data.settings.sig_x*data.settings.sig_y)^(1/2);
data.tau = getTauExp(data.sig0,data.Ti+data.Te);

% determine ion variable names based on equation set
if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
    data.i.n = 'n';
    data.i.v_x = 'v_x';
    data.i.v_y = 'v_y';
    data.i.T = 'temp';
    data.e.n = 'n';
    data.e.v_x = 'v_x';
    data.e.v_y = 'v_y';
    data.e.T = 'temp';
elseif max(strcmp({'ideal_2F'},data.config.eq_set))
    data.i.n = 'i_n';
    data.i.v_x = 'i_v_x';
    data.i.v_y = 'i_v_y';
    data.i.T = 'i_temp';
    data.e.n = 'e_n';
    data.e.v_x = 'e_v_x';
    data.e.v_y = 'e_v_y';
    data.e.T = 'e_temp';
elseif strcmp(data.config.eq_set,{'ideal_mhd_2E'})
    data.i.n = 'n';
    data.i.v_x = 'v_x';
    data.i.v_y = 'v_y';
    data.i.T = 'i_temp';
    data.e.n = 'n';
    data.e.v_x = 'v_x';
    data.e.v_y = 'v_y';
    data.e.T = 'e_temp';
end

end