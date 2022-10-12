function [data] = loadData(folder,flags)
% load various files
f = filesep;
data.folder = folder;
data.settings = readSettingsFile([folder f 'plasma.settings']);
data.config = readConfigFile([folder f 'ucnp.config']);
data.state.init = readStateFile([folder f 'init.state'],flags.ghost_cells,data.config.eq_set);
data.grids = readOutFile([folder f 'mhd.out'],flags.ghost_cells);

% determine central electron and ion temperatures from settings file anc oompute tau_exp4
data.Te = data.settings.Te;
if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
    data.Ti = data.settings.Te;
elseif max(strcmp({'ideal_2F','ideal_mhd_2E'},data.config.eq_set))
    data.Ti = data.settings.Ti;
else
    error('Error: equation set is not valid.')
end
data.sig0 = (data.settings.sig_x*data.settings.sig_y^2)^(1/3);
data.tau = getTauExp(data.sig0,data.Ti+data.Te);

end