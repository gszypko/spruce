function [mat_smoothed] = sgfilt2D(mat,kx,ky,px,py)
mat_smoothed = mat;
bx = floor(kx/2);
by = floor(ky/2);
coeffs = sgsf_2d(-bx:bx,-by:by,px,py,1);

for i = 1+bx:size(mat,1)-bx
    for j = 1+by:size(mat,2)-by
        submat = mat(i-by:i+by,j-bx:j+bx);
        mat_smoothed(i,j) = sum(submat.*coeffs,'all');
    end
end

end