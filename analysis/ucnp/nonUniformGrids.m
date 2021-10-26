clc, clearvars -except inp, close all, f = filesep;

sigx = 0.1;
xmax = 50*sigx;
Nx = 151;
dxmin = sigx/20;

f = @(x) (x-1)^10;

Nxbar = 0;
for i = 1:Nx
    Nxbar = Nxbar + f(i);
end
m = (xmax - Nx*dxmin)/Nxbar;

dx = dxmin;
x = 0;
for i = 2:Nx
    dx(i) = dxmin+m*f(i);
    x(i) = x(i-1)+dx(i-1)/2+dx(i)/2;
end

fig = figure;
dx = [dx(2:end) dx];
x = [-x(2:end) x];
plot(x,dx,'.')
hold on
plot([min(x) max(x)],[sigx/10 sigx/10])
plot(3.*[sigx sigx],[min(dx) max(dx)])

%% use exponential growth