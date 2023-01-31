clc, clear, setpath; close all

N = 301; % total number of grid cells along a given dimension
sig = 0.1; % RMS plasma size
r_max = 6*sig; % the spatial extent of the grid is [-1 1]*r_max
A = 1.07; % multiplicative growth factor for grid cells
B = 0.99999; % secondary growth factor
[x,dx] = getNonUniformGrids(N,r_max,A,B,true);
x = x./sig;
dx = dx./sig;
fig = figure;
fig.Color = [1 1 1];
fig.Position = [612   441   430   332];

yyaxis left
plot(x,dx,'.')
hold on
y = [0.1 0.1];
plot([min(x) max(x)],[0.1 0.1])
plot([3 3 -3 -3],[0 1 1 0])

xlabel('r / \sigma')
ylabel('dr / \sigma')
ylim([0 max([dx y])*1.15])
grid minor

yyaxis right
plot(x(2:end),abs(diff(dx)./dx(2:end)),'.')
xlim([min(x) max(x)])
ylim([0 .1])
ylabel('Fractional Change in dr')

dlm = ' - ';
str1 = ['A = ' num2str(A)];
str2 = ['B = ' num2str(B)];
title([str1 dlm str2],'FontWeight','normal')
