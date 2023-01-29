function [data] = loadData(folder,flags)
% load settings, config, and state files
f = filesep;
data.folder = folder;
data.settings = readSettingsFile([folder f 'plasma.settings']);
data.config = readConfigFile([folder f 'ucnp.config']);
data.state.init = readStateFile([folder f 'init.state'],flags.ghost_cells,data.config.eq_set);

% determine grids to be loaded from mhd.out
grids_to_load = {''};
if flags.plot_grids
    grids_to_load = [grids_to_load flags.vars]; 
end
if flags.vlasov_analysis
    if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
        vlasov_grids = {'n','v_x','v_y','temp'};
    elseif max(strcmp({'ideal_2F','ideal_mhd_2E'},data.config.eq_set))
        vlasov_grids = {'i_n','i_v_x','i_v_y','i_temp','e_temp'};
    elseif strcmp(data.config.eq_set,{'ideal_mhd_2E'})
        vlasov_grids = {'n','v_x','v_y','i_temp','e_temp'};
    end
    grids_to_load = [grids_to_load vlasov_grids];
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
data.sig0 = (data.settings.sig_x*data.settings.sig_y)^(1/2);
data.tau = getTauExp(data.sig0,data.Ti+data.Te);

end