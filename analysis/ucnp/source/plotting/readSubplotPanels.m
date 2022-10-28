function [panel] = readSubplotPanels()
gitdir = setpath(false);
f = filesep;
[C] = readfile([gitdir f 'source' f 'plotting' f 'subplot_panels.txt']);
Csplit = strtrim(split(C,'='));
panel.num = str2double(Csplit(:,1));
panel.row = str2double(extractBefore(Csplit(:,2),'x'));
panel.col = str2double(extractAfter(Csplit(:,2),'x'));
panel.pos = str2double(split(Csplit(:,3),','));
end