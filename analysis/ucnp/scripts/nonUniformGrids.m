clc, clear, close all

N = 301; % total number of grid cells along a given dimension
sig = 0.1; % RMS plasma size
r_max = 7.5*sig; % the spatial extent of the grid is [-1 1]*r_max
A = 1.06; % multiplicative growth factor for grid cells
B = .99994; % secondary growth factor
dx = 1; % start value for growth equation
x = 0;
for i = 2:ceil(N/2)
    dx(i) = 1 + A*(dx(i-1)-B);
    x(i) = x(i-1)+dx(i-1)/2+dx(i)/2;
end
dx = [dx(end:-1:2) dx]./max(x)*r_max./sig;
x = [-x(end:-1:2) x]./max(x)*r_max./sig;

fig = figure;
fig.Color = [1 1 1];
fig.Position = [612   441   430   332];

yyaxis left
plot(x,dx,'.')
hold on
y = [0.1 0.1];
plot([min(x) max(x)],[0.1 0.1])
plot([3 3],[0 1])

xlabel('r / \sigma')
ylabel('dr / \sigma')
ylim([0 max([dx y])*1.15])
grid minor

yyaxis right
plot(x(2:end),abs(diff(dx)./dx(2:end)),'.')
ylim([0 .1])
ylabel('Fractional Change in dr')

dlm = ' - ';
str1 = ['A = ' num2str(A)];
str2 = ['B = ' num2str(B)];
title([str1 dlm str2],'FontWeight','normal')
saveas(fig,'fig.png')