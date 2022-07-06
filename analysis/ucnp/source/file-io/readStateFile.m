function [out] = readStateFile(path,removeGhostCells,Ng,eq_set)
% directory (string): full path to directory containing 'plasma.settings'
% opt (bool): (true) trim ghost cells from matrices (false) do not
% out (struct): fields of <out> contain plasma quantities in cgs units, see below

if ~removeGhostCells, Ng = 0; end

% ensure that .state file exists
if ~endsWith(path,'.state'), error('File extension must be .state'); end

C = readfile(path);
dlm = ',';

temp = split(C{2},dlm);
out.Nx = str2double(temp{1}) - 2*Ng;
out.Ny = str2double(temp{2}) - 2*Ng;
out.mI = str2double(C{4});
out.index = str2double(C{6});
if strcmp(eq_set,'ideal_mhd')
    q = {'d_x','d_y','pos_x','pos_y','rho','temp','mom_x','mom_y','be_x','be_y'};
elseif strcmp(eq_set,'ideal_mhd_2E')
    q = {'d_x','d_y','pos_x','pos_y','rho','i_temp','e_temp','mom_x','mom_y','be_x','be_y'};
else
    error('Error: equation set is not valid.')
end
qind = zeros(size(q));
for i = 1:length(C)
    ind = find(strcmp(C{i},q));
    if ~isempty(ind)
        qind(ind) = i;
    end
end

for i = 1:length(q)
    out.grids.(q{i}) = getgrid(C,qind(i),Ng,out.Nx,out.Ny);
end

out.pos_x = out.grids.pos_x(1,:);
out.pos_y = out.grids.pos_y(:,1)';

end