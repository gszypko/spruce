clearvars
datadir = 'C:\Users\grant\Documents\GitHub\mhd\output';
f = filesep;
folders = dir(datadir);
folders(1:2) = [];
for i = 1:length(folders)
    tok_before = extractBefore(folders(i).name,'_');
    tok_after = extractAfter(folders(i).name,'_');
    num = str2double(tok_after);
    num = num + 27;
    new_name = [tok_before '_' num2str(num)];
    movefile([folders(i).folder f folders(i).name],[folders(i).folder f new_name]);
end