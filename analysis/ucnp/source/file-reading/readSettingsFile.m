function [out] = readSettingsFile(directory)
% directory (string): full path to directory containing 'plasma.settings'
% out (struct): fields of <out> contain plasma quantities in cgs units, see below

% ensure that plasma.settings file exists
f = filesep;
files = dir(directory);
ind = strcmp({files.name},'plasma.settings');
if max(ind)~=1, error('The file <plasma.settings> not found in <directory>.'); end

% read in contents from file
fid = fopen([files(ind).folder f files(ind).name],'r');
iter = 0;
while ~feof(fid)
    iter = iter + 1;
    line = fgetl(fid);
    C{iter} = strsplit(line,',');
end

names = C{1};
units = C{2};
vals = C{3};

% reformat file contents into output structure
out = struct;
q = {'n','sig_x','sig_y','Ti','Te','dBdx','x_lim','y_lim','Nx','Ny','m_i','gam'};
qstr = {'n','sigx','sigy','Ti','Te','dBdx','xlim','ylim','Nx','Ny','mI','adiabatic_index'};
for i = 1:length(names)
    ind = find(strcmp(names{i},q)); 
    if length(ind) == 1
        out.(qstr{ind}) = str2double(vals{i});
        
        ind2 = find(strcmp(units{i},names));
        if length(ind2) == 1
            if ~strcmp(units{ind2},'cgs')
                error('You can only express quantities in terms of other quantities that are in cgs units.')
            end
            out.(qstr{ind}) = out.(qstr{ind})*str2double(vals{ind2});
        end
    end

end

end