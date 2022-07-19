%% Initiate Program
clc, clearvars -except inp, close all, f = filesep; setpath;

if ~exist('inp','var')
    %% Generate Input Cell
    inp = {'folder'};
    openvar('inp')
end
fields = {'folder'};
s = spreadsheet2struct(inp,fields);

% user controls
plotGridTimeEvol = true;
doGaussianAnalysis = true;
removeGhostCells = true;
numGhostCells = 2;
eic_opt = true;

%% Read in and Analyze Data
% define constnats
c = defineConstants();

% read .settings file
disp('Starting Analysis...')
for i = 1:length(s)
    disp(['Data set: ' num2str(i) '/' num2str(length(s))])
    data = loadData(s(i).folder,removeGhostCells,numGhostCells);
    if plotGridTimeEvol, plotGridEvol(data); end
    if doGaussianAnalysis, gaussianAnalysis(data,eic_opt); end
end
disp('Analysis Complete.')