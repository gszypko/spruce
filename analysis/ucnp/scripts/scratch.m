clc, clearvars, f = filesep;
maindir = 'C:\Users\Grant\OneDrive\Research\mhd\data-expsims\an-mhd-09.26.22';
s = dir(maindir);
s(1:2) = [];
for i = 1:length(s)
    s(i).s = dir([s(i).folder f s(i).name]);
    s(i).s(1:2) = [];
    for j = 1:length(s(i).s)
        name_to_delete = [s(i).s(j).folder f s(i).s(j).name];
        isdir = s(i).s(j).isdir;
        if ~endsWith(name_to_delete,'os.mat')
            disp(name_to_delete);
            if ~isdir
                delete(name_to_delete);
            else
                rmdir(name_to_delete,'s');
            end
        end
    end
end