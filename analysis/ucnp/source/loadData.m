function [data] = loadData(folder,flags)
% load settings, config, and state files
f = filesep;
data.folder = folder;
data.settings = readSettingsFile([folder f 'plasma.settings']);
data.config = readConfigFile([folder f 'ucnp.config']);
data.state.init = readStateFile([folder f 'init.state'],flags.ghost_cells,data.config.eq_set);

% determine adiabatic indices for each species
data.adiabatic_index_e = data.settings.adiabatic_index;
data.adiabatic_index_i = data.settings.adiabatic_index;
if data.config.global_temp_e, data.adiabatic_index_e = 1; end
if data.config.global_temp_i, data.adiabatic_index_i = 1; end

% determine grids to be loaded from mhd.out
grids_to_load = {''};
if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
    base_grids = {'n','v_x','v_y','temp'};
elseif max(strcmp({'ideal_2F'},data.config.eq_set))
    base_grids = {'i_n','i_v_x','i_v_y','i_temp','e_temp'};
elseif strcmp(data.config.eq_set,{'ideal_mhd_2E'})
    base_grids = {'n','v_x','v_y','i_temp','e_temp'};
end
grids_to_load = [grids_to_load base_grids];
if flags.plot_grids
    grids_to_load = [grids_to_load flags.vars]; 
end

grids_to_load = grids_to_load(~strcmp(grids_to_load,''));
grids_to_load = unique(grids_to_load);

% load grids from mhd.out
data.grids = readOutFile([folder f 'mhd.out'],flags.ghost_cells,grids_to_load,flags.time_window,flags.time_interval);

% determine central electron and ion temperatures from settings file anc compute tau_exp
data.Te = data.settings.Te;
if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
    data.Ti = data.settings.Te;
elseif max(strcmp({'ideal_2F','ideal_mhd_2E'},data.config.eq_set))
    data.Ti = data.settings.Ti;
else
    error('Error: equation set is not valid.')
end
data.sigx = data.settings.sig_x;
data.sigy = data.settings.sig_y;
data.sig0 = (data.settings.sig_x*data.settings.sig_y)^(1/2);
data.tau = getTauExp(data.sig0,data.Ti+data.Te);

end