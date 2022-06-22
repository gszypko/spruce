%% Initiate Program
clc, clearvars -except inp, close all, f = filesep; setpath;

if ~exist('inp','var')
    %% Generate Input Cell
    inp = {};
    openvar('inp')
end
fields = {'date','phase','fields','Te','n_dist','n_max','n_min','sigx','sigy','grid_opt','xlim','ylim','tmax','Nx',...
    'Ny','bound_strength','bound_decay','visc_strength','set','folder'};
if size(inp,2) ~= length(fields), error('The number of columns in ''inp'' must match the length of ''fields''.'); end
s = cell2struct(inp,fields,2);

% user controls
plotGridTimeEvol = true;
doGaussianAnalysis = true;
removeGhostCells = true;
numGhostCells = 2;

%% Read in and Analyze Data
% define constnats
c = defineConstants();

% read .settings file
disp('Starting Analysis...')
for i = 1:length(s)
    disp(['Data set: ' num2str(i) '/' num2str(length(s))])
    data = loadData(s(i).folder,removeGhostCells,numGhostCells);
    if plotGridTimeEvol, plotGridEvol(data); end
    if doGaussianAnalysis, gaussianAnalysis(data); end
end
disp('Analysis Complete.')