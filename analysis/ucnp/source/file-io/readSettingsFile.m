function [out] = readSettingsFile(path)
% directory (string): full path to directory containing 'plasma.settings'
% out (struct): fields of <out> contain plasma quantities in cgs units, see below

% ensure that plasma.settings file exists
if ~endsWith(path,'.settings'), error('<path> does not point to .settings file'); end

% read in contents from file
C = readfile(path);

names = cell(length(C),1);
vals = cell(length(C),1);
units = cell(length(C),1);
for i = 1:length(C)
    line = split(C{i},"=");
    names{i} = strtrim(line{1});
    units{i} = strtrim(line{2});
    vals{i} = strtrim(line{3});
end

% reformat file contents into output structure
out = struct;
q = {'n','sig_x','sig_y','Ti','Te','dBdx','x_lim','y_lim','Nx','Ny','m_i','gam'};
qstr = {'n','sigx','sigy','Ti','Te','dBdx','xlim','ylim','Nx','Ny','mI','adiabatic_index'};
specified = zeros(size(q));
for i = 1:length(names)
    ind = find(strcmp(names{i},q)); 
    if length(ind) == 1
        specified(ind) = 1;
        if strcmp(units{i},'opt')
            out.(qstr{ind}) = vals{i};
        else
            out.(qstr{ind}) = str2double(vals{i});
        end
        
        if ~max(strcmp(units{i},{'cgs','opt'}))
            ind2 = find(strcmp(units{i},names));
            if isempty(ind2)
                error('Variable unit must be ''cgs'', ''opt'', or another variable name.')
            elseif length(ind2) == 1
                if ~strcmp(units{ind2},'cgs')
                    error('You can only express quantities in terms of other quantities that are in cgs units.')
                end
                out.(qstr{ind}) = out.(qstr{ind})*str2double(vals{ind2});
            end
        end
    end

end

end