function [data] = loadData(folder,removeGhostCells,numGhostCells)

f = filesep;

data.folder = folder;
data.settings = readSettingsFile([folder f 'plasma.settings']);
data.state.init = readStateFile([folder f 'init.state'],removeGhostCells,numGhostCells);
data.state.end = readStateFile([folder f 'end.state'],removeGhostCells,numGhostCells);
data.grids = readOutFile([folder f 'mhd.out'],removeGhostCells,numGhostCells);

fields = fieldnames(data.grids.vars);
if max(strcmp(fields,'rho'))
    for i = 1:length([data.grids.vars.time])
        data.grids.vars(i).n = data.grids.vars(i).rho./data.settings.mI;
    end
end

Te = data.settings.Te;
data.sig0 = (data.settings.sigx*data.settings.sigy^2)^(1/3);
data.tau = getTauExp(data.sig0,Te+Te);

end