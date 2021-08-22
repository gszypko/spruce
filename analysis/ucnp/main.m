%% Initiate Program
clc, clear, close all, f = filesep;

% include MHD analysis files for UCNPs in Matlab search path
andir = extractBefore(matlab.desktop.editor.getActiveFilename,mfilename);
gitdir = extractBefore(andir,'analysis/ucnp');
addpath(genpath(andir))

% full path to MHD simulation files
datadir = [gitdir f 'output' f 'array0'];
figdir = [andir f 'figs'];
mkdir(figdir)

%% 

data = readStateFile(datadir);

varnames = data.grid_names;

for i = 1:length(varnames)    
    fig = figure;
    fig.Position = [490   283   540   420];
    fig.Color = [1 1 1];

    ax = axes();

    xdata = [data.x_vec];
    ydata = [data.y_vec];
    zdata = [data.(varnames{i})];
    imagesc(xdata,ydata,zdata)

    cb = colorbar;
    cb.FontSize = 12;
    cb.Label.String = varnames{i};
    ax.YDir = 'Normal';

    ax.PlotBoxAspectRatio = [1 1 1];
    ax.FontSize = 12;
    ax.XTick = round(linspace(min(xdata),max(xdata),5),2,'significant');
    ax.YTick = round(linspace(min(ydata),max(ydata),5),2,'significant');
    
    xlabel('x (cm)')
    ylabel('y (cm)')
    
    savepath = [figdir f 'var' num2str(i) '-' varnames{i} '.png'];
    saveas(fig,savepath)
    close(fig)
end