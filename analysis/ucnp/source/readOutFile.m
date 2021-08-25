function [out] = readOutFile(directory,removeGhostCells)
% directory (string): full path to folder containing MHD simulation files (i.e., .state, .settings, .config)
% removeGhostCells (boolean): (true) trim ghost cells from matrix (false) do not

% Note: all units in the .state file are in cgs

%% Read MHD .state File
% define # of ghost cells
numGhostCells = 2;

% specify full path to .state file
f = filesep;
filename = 'mhd.out';
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

% begin parsing non-grid quantities
temp = split(data{2},dlm);
out.xdim = str2double(temp{1}); % number of elements in each grid col
out.ydim = str2double(temp{2}); % number of elements in each grid row
out.pos_x = gridvec2mat(data{4},removeGhostCells,numGhostCells);
out.pos_y = gridvec2mat(data{6},removeGhostCells,numGhostCells);
out.b_x = gridvec2mat(data{8},removeGhostCells,numGhostCells);
out.b_y = gridvec2mat(data{10},removeGhostCells,numGhostCells);
out.b_z = gridvec2mat(data{12},removeGhostCells,numGhostCells);

% get time values
ind = find(startsWith(data,'t='))';
time = zeros(size(ind));
rowsPerT = ind(2) - ind(1);
out.vars = struct;
out.vars(length(time)).time = [];
for i = 1:length(ind)
    out.vars(i).time = str2double(extractAfter(data{ind(i)},'t='));
    for j = (ind(i)+1):2:(ind(i)+rowsPerT-2)
        out.vars(i).(data{j}) = gridvec2mat(data{j+1},removeGhostCells,numGhostCells);
    end
end
out.xVec = out.pos_x(1,:);
out.yVec = out.pos_y(:,1)';

gridNames = fieldnames(out.vars);
out.gridNames = gridNames(2:end);

end