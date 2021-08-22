function [out] = readStateFile(directory)
% directory (string): full path to folder containing MHD simulation files (i.e., .state, .settings, .config)

% Note: all units in the .state file are in cgs

%% Read MHD .state File
% specify full path to .state file
f = filesep;
filename = 'mhd.state';
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
dlm1 = ','; % element delimiter
dlm2 = ';'; % row delimiter for matrix

% begin parsing non-grid quantities
temp = split(data{2},dlm1);
out.xdim = str2double(temp{1}); % number of elements in each grid col
out.ydim = str2double(temp{2}); % number of elements in each grid row
out.ion_mass = str2double(data{4}); % ion mass (g)
out.adiabatic_index = str2double(data{6}); % adiabatic index for ion

% begin parsing grid quantities
    % In the MHD code, a matrix M(i,j) is set up such that i == row == x axis and j == column == y axis.
    % This means x position varies with row number and vice versa.
    % For visualization using 'imagesc', I want x position to correspond to column number.
    % In the following parsing, I reverse the row/column association with x/y axis that is contained within the .state file.
grid_names = {'pos_x','pos_y','rho','temp','mom_x','mom_y','b_x','b_y','b_z'};
out.grid_names = grid_names;
for i = 1:length(grid_names)
    ind = find(strcmp(data,grid_names{i}));
    out.(grid_names{i}) = zeros(out.xdim,out.ydim);
    temp_row = split(data{ind+1},dlm2,2); % each element of 'temp_row' contains a row vector
    for j = 1:length(temp_row)
        temp = split(temp_row{j},dlm1,2);
        out.(grid_names{i})(j,:) = str2double(temp);
    end
    out.(grid_names{i}) = out.(grid_names{i})';
end
out.x_vec = out.pos_x(1,:);
out.y_vec = out.pos_y(:,1)';

% trime off ghost cells - first two and last two rows and columns are removed
out.xdim = out.xdim - 2;
out.ydim = out.ydim - 2;
for i = 1:length(grid_names)
    out.(grid_names{i}) = out.(grid_names{i})(3:end-2,3:end-2);
end

