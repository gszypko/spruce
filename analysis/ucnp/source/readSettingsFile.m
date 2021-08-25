function [out] = readSettingsFile(directory)
% directory (string): full path to MHD simulation folder, containing .settings file

% define path to .settings file
f = filesep; 
filepath = [directory f 'plasma.settings']; 

% read in each line from .settings file
fileID = fopen(filepath);
iter = 0;
data = cell(1000,1);
while ~logical(feof(fileID))
    iter = iter + 1;
    data{iter,1} = fgetl(fileID);
end
data = data(1:iter);

% parse comma delimited lines
dlm = ',';
names = split(data(1),dlm);
units = split(data(2),dlm);
values = split(data(3),dlm);

out = struct;
out(length(names)).name = [];
for i = 1:length(names)
    out(i).name = names{i};
    out(i).unit = units{i};
    out(i).valstr = values{i};
    out(i).valnum = str2double(values{i});
end

end

