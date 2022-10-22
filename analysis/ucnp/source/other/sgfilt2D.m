function [mat_sg,ind_x,ind_y] = sgfilt2D(mat,kx,ky,px,py,opt)
% mat (matrix): matrix to be smoothed
% kx (double): kernel size along second (x) dimension - columns
% ky (double): kernel size along first (y) dimension - rows
% px (double): polynomial order for second dimension
% py (double): polynomial order for first dimension
% opt (boolean): (true) remove boundary cells (false) leave boundary cells
% ind_x (vector): indices corresponding to x interior cells
% ind_y (vector): indices corresponding to y interior cells
% mat_sg (matrix): smoothed matrix, with or without boundary cells depending on <opt>

if nargin < 6, opt = false; end

% determine number of boundary cells and coefficient matrix
bx = floor(kx/2);
by = floor(ky/2);
coeffs = sgsf_2d(-bx:bx,-by:by,px,py,1);

% compute sg smoothed value for each interior cell
mat_sg = mat;
ind_x = 1+bx:size(mat,2)-bx;
ind_y = 1+by:size(mat,1)-by;    
for i = ind_y
    for j = ind_x
        submat = mat(i-by:i+by,j-bx:j+bx);
        mat_sg(i,j) = sum(submat.*coeffs,'all');
    end
end

% option for removing boundary cells from matrix
if opt
    mat_sg = mat_sg(1+by:end-by,1+bx:end-bx);
end

end