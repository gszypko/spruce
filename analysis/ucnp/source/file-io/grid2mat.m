function [mat] = grid2mat(grid,Nx,Ny,opt)

mat = zeros(Ny,Nx);
iter = 0;
for i = 1:Nx
    for j = 1:Ny
        iter = iter + 1;
        mat(j,i) = str2double(grid(iter));
    end
end

if opt
    mat = mat(3:end-2,3:end-2);
end

end