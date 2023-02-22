clc, clear, close all

n_max = 1e9;
n_min = 1e4;
N = 301;
sig = 0.1;
lim = 10*sig;

[x,dx] = getNonUniformGrids(N,lim,1.05,.9999,true);
[y,dy] = getNonUniformGrids(N,lim,1.05,.9999,true);
[X,Y] = meshgrid(x,y);
n = n_max*exp(-abs(X)./sig) + n_min;
dndx = derivative1D(dx,dy,n,'x');
dndy = derivative1D(dx,dy,n,'y');

[~,center_ind] = min(abs(y));
plot(x,dndx(center_ind,:)./n(center_ind,:),'.');
yyaxis right
plot(x,n(center_ind,:),'.');

s.n = n;
s.dndx = dndx;
s.dndy = dndy;

vars = {'n','dndx','dndy'};
varstr = {'n','dndx','dndy'};
num = length(vars);
[fig,ax,an,row,col] = open_subplot(num);

iter = 0;
for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        iter = iter + 1;
        if iter > num, break, end

        cax = get_axis(fig,ax{i,j});
        xdata = x;
        ydata = y;
        zdata = s.(vars{iter})./s.n;
        imagesc(xdata,ydata,zdata)
        colorbar
        cax.YDir = 'Normal';
        cax.PlotBoxAspectRatio = [1 1 1];
        cax.FontSize = 10;
        if i == size(ax,1), xlabel('x (cm)'), end
        if j == 1, ylabel('y (cm)'), end
        title(varstr{iter},'FontWeight','normal')
    end
end
