function [os] = loadData(folder,removeGhostCells)

os.folder = folder;
os.settings = readSettingsFile(folder);
os.state = readStateFile(folder,removeGhostCells);
os.grids = readOutFile(folder,removeGhostCells);

fields = fieldnames(os.grids.vars);
if max(strcmp(fields,'rho'))
    for i = 1:length(os.grids.time)
        os.grids.vars(i).n = os.grids.vars(i).rho./os.settings.mI;
    end
end

end