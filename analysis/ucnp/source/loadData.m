function [data] = loadData(folder,removeGhostCells,numGhostCells)

f = filesep;

data.folder = folder;
data.settings = readSettingsFile([folder f 'plasma.settings']);
data.config = readConfigFile([folder f 'ucnp.config']);
data.state.init = readStateFile([folder f 'init.state'],removeGhostCells,numGhostCells,data.config.eq_set);
data.state.end = readStateFile([folder f 'end.state'],removeGhostCells,numGhostCells,data.config.eq_set);
data.grids = readOutFile([folder f 'mhd.out'],removeGhostCells,numGhostCells);

fields = fieldnames(data.grids.vars);
if max(strcmp(fields,'rho'))
    for i = 1:length([data.grids.vars.time])
        data.grids.vars(i).n = data.grids.vars(i).rho./data.settings.mI;
    end
end

data.Te = data.settings.Te;
if strcmp(data.config.eq_set,'ideal_mhd') || strcmp(data.config.eq_set,'ideal_mhd_cons')
    data.Ti = data.settings.Te;
elseif strcmp(data.config.eq_set,'ideal_mhd_2E')
    data.Ti = data.settings.Ti;
else
    error('Error: equation set is not valid.')
end
data.sig0 = (data.settings.sigx*data.settings.sigy^2)^(1/3);
data.tau = getTauExp(data.sig0,data.Ti+data.Te);

end