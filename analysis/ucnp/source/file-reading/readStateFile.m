function [out] = readStateFile(directory,opt)
% directory (string): full path to directory containing 'plasma.settings'
% opt (bool): (true) trim ghost cells from matrices (false) do not
% out (struct): fields of <out> contain plasma quantities in cgs units, see below

% ensure that plasma.settings file exists
f = filesep;
files = dir(directory);
ind = strcmp({files.name},'init.state');
if max(ind)~=1, error('The file <mhd.state> not found in <directory>.'); end

% read in contents from file
fid = fopen([files(ind).folder f files(ind).name],'r');
iter = 0;
while ~feof(fid)
    iter = iter + 1;
    line = fgetl(fid);
    C{iter} = strsplit(line,{',',';'});
end

out.Nx = str2double(C{2}{1});
out.Ny = str2double(C{2}{2});
out.mI = str2double(C{4}{1});
out.index = str2double(C{6}{1});

q = {'d_x','d_y','pos_x','pos_y','rho','temp','mom_x','mom_y','b_x','b_y'};
qind = [];
for i = 1:length(C)
    ind = find(strcmp(C{i}{1},q));
    if ~isempty(ind)
        qind(ind) = i+1;
    end
end

for i = 1:length(q)
    out.grids.(q{i}) = grid2mat(C{qind(i)},out.Nx,out.Ny,opt);
end

out.pos_x = out.grids.pos_x(1,:);
out.pos_y = out.grids.pos_y(:,1)';

end