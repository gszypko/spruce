% transect for plot
x = linspace(-.5,.5,301);
y = 0;

% grid data
X = data.grids.pos_x;
Y = data.grids.pos_y;
n = data.grids.vars(10).n;
n = sgolayfilt(n,3,25);

% interpolation
n_int = interp2(X,Y,n,x,y);
n_grad = diff(n_int);
% n_grad = sgolayfilt(n_grad,3,25);
n_int_less = zeros(1,length(n_grad));
for i = 1:length(n_int_less)
    n_int_less(i) = mean(n_int(i:i+1));
    x_less(i) = mean(x(i:i+1));
end

gradnOn = n_grad./n_int_less;

figure
plot(x_less(25:end-25),gradnOn(25:end-25))
