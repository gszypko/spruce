function [data] = loadData(folder,numGhostCells)

f = filesep;

data.folder = folder;
data.settings = readSettingsFile([folder f 'plasma.settings']);
data.config = readConfigFile([folder f 'ucnp.config']);
data.state.init = readStateFile([folder f 'init.state'],numGhostCells,data.config.eq_set);
% data.state.end = readStateFile([folder f 'end.state'],numGhostCells,data.config.eq_set);
data.grids = readOutFile([folder f 'mhd.out'],numGhostCells);

fields = fieldnames(data.grids.vars);
if max(strcmp(fields,'rho'))
    for i = 1:length([data.grids.vars.time])
        data.grids.vars(i).n = data.grids.vars(i).rho./data.settings.mI;
    end
end

data.Te = data.settings.Te;
if max(strcmp({'ideal_mhd','ideal_mhd_cons'},data.config.eq_set))
    data.Ti = data.settings.Te;
elseif max(strcmp({'ideal_2F','ideal_mhd_2E'},data.config.eq_set))
    data.Ti = data.settings.Ti;
else
    error('Error: equation set is not valid.')
end
data.sig0 = (data.settings.sigx*data.settings.sigy^2)^(1/3);
data.tau = getTauExp(data.sig0,data.Ti+data.Te);

end