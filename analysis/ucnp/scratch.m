clc, clearvars -except folders
f = filesep;
newfolders = cell(size(folders));

%% update folder names
for i = 1:length(folders)
    currname = folders{i};
    % remove apostrophe from name
    if startsWith(currname,'''')
        currname = currname{2:end-1};
    end
    relfolder = extractAfter(currname,'C:\');
    relfolder = strrep(relfolder,'/',f);
    newname = ['C:\Users\grant\OneDrive\Research\mhd' f relfolder];
    if endsWith(newname,f)
        newname(end) = [];
    end
    newname = strrep(newname,'\\',f);
    disp(newname)
end

%% remove unnecessary analysis subfolders

% loop through data folders
for i = 1:length(folders)
    
    % check that folder exists
    if ~exist(folders{i}, 'dir')
        disp(folders{i})
    end
    
    % loop through subfolders
    subfolders = dir(folders{i});
    subfolders = subfolders(3:end);
    subfolders(~[subfolders.isdir]) = [];
    if isempty(subfolders)
        continue;
    end
    for j = 1:length(subfolders)
        if ~strcmp(subfolders(j).name,'an-RELACS-10.13.21')
            dirtodelete = [subfolders(j).folder f subfolders(j).name];
            rmdir(dirtodelete,'s')
        end
    end
end