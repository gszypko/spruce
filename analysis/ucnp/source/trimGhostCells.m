function [matOut] = trimGhostCells(matIn)
% matIn (m x n double): matrix imported from MHD simulation, including ghost cells
% matOut (m-2 x n-2 double): matrix with ghost cells removed

% This function removes the first two rows and columns from 'matIn'

matOut = matIn(3:end-2,3:end-2);
end

