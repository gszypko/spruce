function [x,y,mat_sg] = sgfilt2D(x,y,mat,kx,ky,px,py,opt)
% mat (mxn double)
% kx (double): kernel size along second (x) dimension - columns
% ky (double): kernel size along first (y) dimension - rows
% px (double): polynomial order for second dimension
% py (double): polynomial order for first dimension
% opt (boolean): (true) remove boundary cells (false) leave boundary cells

% determine number of boundary cells and coefficient matrix
bx = floor(kx/2);
by = floor(ky/2);
coeffs = sgsf_2d(-bx:bx,-by:by,px,py,1);

% compute sg smoothed value for each interior cell
mat_sg = mat;
for i = 1+bx:size(mat,1)-bx
    for j = 1+by:size(mat,2)-by
        submat = mat(i-by:i+by,j-bx:j+bx);
        mat_sg(i,j) = sum(submat.*coeffs,'all');
    end
end

% option for removing boundary cells from matrix
if opt
    x = x(1+bx:end-bx);
    y = y(1+by:end-by);
    mat_sg = mat_sg(1+by:end-by,1+bx:end-bx);
end

end