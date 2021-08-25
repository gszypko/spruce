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
mkdir([figdir f 'state'])
mkdir([figdir f 'out'])

%% Read in Data
% define constnats
c = constants();

% read .settings file
plasma = readSettingsFile(datadir);

% read MHD simulation data
removeGhostCells = true;
fileType = 'out';
if strcmp(fileType,'state')
    mhd = readStateFile(datadir,removeGhostCells);
elseif strcmp(fileType,'out')
    mhd = readOutFile(datadir,removeGhostCells);
end

%% Plot State Variables
if strcmp(fileType,'state')
    varnames = mhd.gridNames;
    for i = 1:length(varnames)    
        fig = figure;
        fig.Position = [490   283   540   420];
        fig.Color = [1 1 1];

        ax = axes();

        xdata = [mhd.xVec];
        ydata = [mhd.yVec];
        zdata = [mhd.(varnames{i})];
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

        savepath = [figdir f 'state' f 'var' num2str(i) '-' varnames{i} '.png'];
        saveas(fig,savepath)
        close(fig)
    end
end

%% Plot Out Variables
if strcmp(fileType,'out')
    gridNames = {'rho','temp','press','v_x','v_y','thermal_energy'};
    for k = 1:length([mhd.vars.time])
        fig = figure;
        fig.Position = [432         207        1039         525];
        fig.Color = [1 1 1];

        colvar = {''};
        rowvar = {''};
        row = 2;
        col = 3;
        num = row*col;

        ax = cell(row,col);
        iter = 0;
        for i = 1:row
            for j = 1:col
                if iter > num - 1, break, end
                iter = iter + 1;
                ax{i,j} = subplot(row,col,iter);

                xdata = [mhd.xVec];
                ydata = [mhd.yVec];
                zdata = [mhd.vars(k).(gridNames{iter})];
                imagesc(xdata,ydata,zdata)

                cb = colorbar;
                cb.FontSize = 12;
                cb.Label.String = gridNames{iter};
                ax{i,j}.YDir = 'Normal';

                ax{i,j}.PlotBoxAspectRatio = [1 1 1];
                ax{i,j}.FontSize = 12;
                
                if i == row, xlabel('x (cm)'), end
                if j == 1, ylabel('y (cm)'), end

            end
        end

        an = annotation('textbox');
        an.Position = [0.0058    0.9219    0.9919    0.0686];
        an.HorizontalAlignment = 'center';
        an.VerticalAlignment = 'middle';
        an.LineStyle = 'none';
        an.FontSize = 12;
        
        dlm = ' - ';
        str1 = ['Iter = ' num2str(k-1)];
        str2 = ['t = ' num2str(mhd.vars(k).time*1e6) '\mus'];
        an.String = [str1 dlm str2];

        savepath = [figdir f 'out' f 'tEvol' num2str(k-1) '.png'];
        saveas(fig,savepath)
        close(fig)
    end
end

