function [grid_names,grid_str] = readGridNames()

gitdir = setpath(false);
C = readfile([gitdir filesep 'gridnames.txt']);
grid_names = strtrim(extractBefore(C,'='));
grid_str = strtrim(extractAfter(C,'='));

end