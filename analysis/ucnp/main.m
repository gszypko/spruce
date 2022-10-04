%% Initiate Program
clc, clearvars -except inp, close all, f = filesep; setpath;

if ~exist('inp','var')
    %% Generate Input Cell
    inp = {'folder';'C:\Users\Grant\OneDrive\Research\mhd\data-expsims\an-mhd-09.26.22\C4.3\set_3'};
end
fields = {'folder'};
s = spreadsheet2struct(inp,fields);

% user controls
plotGridTimeEvol = true;
doGaussianAnalysis = true;
doCompareExpAndSimData = true;
plotFreq = 1;
numGhostCells = 2;
eic_opt = false;


%% Read in and Analyze Data
disp('Starting Analysis...')
for i = 1:length(s)
    disp(['Data set: ' num2str(i-1) '/' num2str(length(s)-1)])
    data = loadData(s(i).folder,numGhostCells);
    if plotGridTimeEvol, plotGridEvol(data,plotFreq); end
    if doGaussianAnalysis, gaussianAnalysis(data,eic_opt,plotFreq); end
    if doCompareExpAndSimData, compareExpAndSimData(data,eic_opt); end
end
disp('Analysis Complete.')