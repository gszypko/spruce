%% Initiate Program
clc, clearvars -except inp, close all, f = filesep; setpath;

if ~exist('inp','var')
    %% Generate Input Cell
    inp = {'folder';'C:\data-mhd\09.06.22\set_10'};
end
fields = {'folder'};
s = spreadsheet2struct(inp,fields);

% user controls
plotGridTimeEvol = true;
plotFreq = 1;
doGaussianAnalysis = true;
numGhostCells = 2;
eic_opt = false;


%% Read in and Analyze Data
disp('Starting Analysis...')
for i = 1:length(s)
    disp(['Data set: ' num2str(i-1) '/' num2str(length(s)-1)])
    data = loadData(s(i).folder,numGhostCells);
    if plotGridTimeEvol, plotGridEvol(data,plotFreq); end
    if doGaussianAnalysis, gaussianAnalysis(data,eic_opt,plotFreq); end
end
disp('Analysis Complete.')