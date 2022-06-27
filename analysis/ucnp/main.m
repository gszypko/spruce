%% Initiate Program
clc, clearvars -except inp, close all, f = filesep; setpath;

if ~exist('inp','var')
    %% Generate Input Cell
    inp = {'27-Jun-2022'	'phase'	'fields'	'Ti'	'Te'	'n_dist'	'n_max'	'n_min'	'sigx'	'sigy'	'grid_opt'	'x_lim'	'y_lim'	't_max'	'Nx'	'Ny'	'set'	'folder'};
    openvar('inp')
end
fields = {'n_max','Ti','Te','folder'};
s = spreadsheet2struct(inp,fields);

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