function [out] = readStateFile(directory,removeGhostCells)
% directory (string): full path to folder containing MHD simulation files (i.e., .state, .settings, .config)
% removeGhostCells (boolean): (true) trim ghost cells from edges of matrices (false) do not

% Note: all units in the .state file are in cgs

%% Read MHD .state File
% define # of ghost cells
numGhostCells = 2;

% specify full path to .state file
f = filesep;
filename = 'init.state';
filepath = [directory f filename];

% start file reading
data = cell(1000,1);
fileID = fopen(filepath);
iter = 0;
while ~logical(feof(fileID))
    iter = iter + 1;
    data{iter,1} = fgetl(fileID);
end
fclose(fileID);
data = data(1:iter);

out = struct;
dlm = ','; % element delimiter
temp = split(data{2},dlm);
out.xdim = str2double(temp{1}); % number of elements in each grid col
out.ydim = str2double(temp{2}); % number of elements in each grid row
out.ion_mass = str2double(data{4}); % ion mass (g)
out.adiabatic_index = str2double(data{6}); % adiabatic index for ion

if removeGhostCells
    out.xdim = out.xdim - numGhostCells;
    out.ydim = out.ydim - numGhostCells;
end

gridNames = {'pos_x','pos_y','rho','temp','mom_x','mom_y','b_x','b_y','b_z'};
out.gridNames = gridNames;
for i = 1:length(gridNames)
    ind = find(strcmp(data,gridNames{i}));
    out.(gridNames{i}) = gridvec2mat(data{ind+1},removeGhostCells,numGhostCells);
end
out.xVec = out.pos_x(1,:);
out.yVec = out.pos_y(:,1)';

end

