function [out] = readOutFile(directory,opt)
% directory (string): full path to directory containing 'plasma.settings'
% opt (bool): (true) trim ghost cells from matrices (false) do not
% out (struct): fields of <out> contain plasma quantities in cgs units, see below

% ensure that plasma.settings file exists
f = filesep;
files = dir(directory);
ind = strcmp({files.name},'mhd.out');
if max(ind)~=1, error('The file <mhd.state> not found in <directory>.'); end

% read in contents from file
disp(['Reading mhd.out: '])
fid = fopen([files(ind).folder f files(ind).name],'r');
iter = 0;
while ~feof(fid)
    iter = iter + 1;
    line = fgetl(fid);
    C{iter} = strsplit(line,{',',';'});
end

out.Nx = str2double(C{2}{1});
out.Ny = str2double(C{2}{2});

out.pos_x = grid2mat(C{4},out.Nx,out.Ny,opt);
out.pos_y = grid2mat(C{6},out.Nx,out.Ny,opt);
out.b_x = grid2mat(C{8},out.Nx,out.Ny,opt);
out.b_y = grid2mat(C{10},out.Nx,out.Ny,opt);

iter_t = 0;
for i = 11:length(C)
    disp(['Line: ' num2str(i-10) '/' num2str(length(C)-10)])
    if startsWith(C{i},'t=')
        iter_t = iter_t + 1;
        out.vars(iter_t).time = str2double(extractAfter(C{i}{1},'t='));
    else
        if length(C{i}) == 1
            out.vars(iter_t).(C{i}{1}) = grid2mat(C{i+1},out.Nx,out.Ny,opt);
        end
    end      
end


end