function [grid_names,grid_str] = readGridNames(path)

gitdir = setpath(false);
f = filesep;
[C] = readfile([gitdir f 'gridnames.txt']);
grid_names = strtrim(extractBefore(C,'='));
grid_str = strtrim(extractAfter(C,'='));

end