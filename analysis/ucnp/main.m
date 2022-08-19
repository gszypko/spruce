%% Initiate Program
clc, clearvars -except inp, close all, f = filesep; setpath;

if ~exist('inp','var')
    %% Generate Input Cell
    inp = {'folder';'C:\Users\Grant\Documents\GitHub\mhd\output\set_0'};
end
fields = {'folder'};
s = spreadsheet2struct(inp,fields);

% user controls
plotGridTimeEvol = true;
doGaussianAnalysis = false;
numGhostCells = 0;
eic_opt = false;

%% Read in and Analyze Data
disp('Starting Analysis...')
for i = 1:length(s)
    disp(['Data set: ' num2str(i) '/' num2str(length(s))])
    data = loadData(s(i).folder,numGhostCells);
    if plotGridTimeEvol, plotGridEvol(data); end
    if doGaussianAnalysis, gaussianAnalysis(data,eic_opt); end
end
disp('Analysis Complete.')