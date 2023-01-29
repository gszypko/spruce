function [grid] = getgrid(C,ind_start,Ng,Nx,Ny)

grid = zeros(Ny,Nx);
grid_begin = ind_start + 1 + Ng;
grid_end = ind_start + Ng + Nx;
iter = 0;
for i = grid_begin:grid_end
    iter = iter + 1;
    temp = split(C{i},',');
    line = sscanf(sprintf(' %s',temp{:}),'%f',[1,Inf]);
    grid(:,iter) = line(Ng+1:end-Ng);
end

end