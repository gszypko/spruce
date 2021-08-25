function [mat] = gridvec2mat(gridvec,removeGhostCells,numGhostCells)
% gridvec (row vector string): vectorized grid from MHD simulation, see grid.hpp
% removeGhostCells (boolean): (true) delete first and last two rows/columns from matrix

%% Convert Vectorized Grid to 2D Matlab Matrix
% delimiters for gridvec
dlm1 = ','; % element delimiter
dlm2 = ';'; % row delimiter for matrix

% parse gridvec into matrix
temp_row = split(gridvec,dlm2,2); % each element of 'temp_row' contains a row vector'
mat = zeros(length(temp_row),length(split(temp_row{1},dlm1,2)));
for j = 1:length(temp_row)
    temp = split(temp_row{j},dlm1,2);
    mat(j,:) = str2double(temp);
end
mat = mat';

if removeGhostCells
    mat = mat(numGhostCells:end-numGhostCells,1+numGhostCells:end-numGhostCells);
end

