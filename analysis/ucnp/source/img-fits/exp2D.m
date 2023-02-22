function [out] = exp2D(c,data)
% c (vector double): fit model coefficients
%   c(1): amplitude
%   c(2): x center
%   c(3): 1/e distance along x axis
%   c(4): y center
%   c(5): 1/e distance along y axis
%   c(6): offset

% data (nx2 double): first column is x positions and second is y positions
% x(i) and y(i) are coordinate pairs to evaluate Gaussian at
x = data(:,1);
y = data(:,2);

% out (column vector double): Gaussian function evaluated at x and y
out = c(1).*exp(-sqrt(((x-c(2))./c(4)).^2 + ((y-c(3))./c(5)).^2)) + c(6);

end