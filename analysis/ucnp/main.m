%% Initiate Program
clc, clearvars -except inp, close all, f = filesep; setpath;

if ~exist('inp','var')
    %% Generate Input Cell
    inp = {};
    openvar('inp')
end
fields = {'date','phase','folder','fields','Te','n_dist','n_max','n_min','sigx','sigy','grid_opt','xlim','ylim','tmax','Nx','Ny','bound_strength','bound_decay','visc_strength'};
if size(inp,2) ~= length(fields), error('The number of columns in ''inp'' must match the length of ''fields''.'); end
s = cell2struct(inp,fields,2);

% user controls
removeGhostCells = true;
loadFromBaseFiles = true;
plotGridTimeEvol = true;
doGaussianAnalysis = false;

%% Read in and Process Data
% define constnats
c = defineConstants();

% read .settings file
disp('Start data proocessing...')
for i = 1:length(s)
    disp(['Data set: ' num2str(i) '/' num2str(length(s))])
    filename = 'os.mat';
    filepath = [s(i).folder f filename];
    
    % check to see whether 'os.mat' was previously saved
    files = dir(s(i).folder);
    found = max(strcmp({files.name},filename));
    
    % if processed data is not found or option
    if loadFromBaseFiles || ~found
        os = loadData(s(i).folder,removeGhostCells);
        save(filepath,'os')
    end
end
disp('Data processing complete.')

%% Run MHD Analysis

for i = 1:length(s)
    disp(['Analyzing set: ' num2str(i) '/' num2str(length(s))])
    load([s(i).folder f 'os.mat'],'os')
    
    if plotGridTimeEvol, plotGridEvol(os); end
    if doGaussianAnalysis, gaussianAnalysis(os); end
end

disp('Done')