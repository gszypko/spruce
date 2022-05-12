function [grid] = getgrid(C,ind_start,Ng,Nx,Ny)

grid = zeros(Ny,Nx);
iter = 0;
for i = ind_start+1+Ng:ind_start+Ng+Nx
    iter = iter + 1;
    temp = split(C{i},',');
    line = sscanf(sprintf(' %s',temp{:}),'%f',[1,Inf]);
    grid(:,iter) = line(Ng+1:end-Ng);
end

end