% add repo to path
maindir = 'C:\Users\grant\Documents\GitHub\mhd\analysis\ucnp'; % where the code sits
addpath(maindir); setpath(maindir); % adds repo full repo to matlab path

% load data
folder = 'C:\Users\grant\OneDrive\Research\projects\23.12.12-UCNP-Shock-Paper-1\Phase-7.11\set_0';
files = dir(folder);
data_mat_found = max(strcmp({files.name},'data.mat'));
if data_mat_found
    load([folder filesep 'data.mat'],'data');
else
    data = loadData(folder);
    save([data.folder filesep 'data.mat'],"data",'-mat');
end

Te = zeros(size(data.grids.time));
for i = 1:length(data.grids.time)
    Te(i) = mean(data.grids.vars(i).e_temp,'all');
end