function [out] = gaussian2D(c,data,offset)
% c (vector double): fit model coefficients
%   c(1): amplitude
%   c(2): x center
%   c(3): rms radius on x axis
%   c(4): y center
%   c(5): rms radius on y axis
%   c(6): offset, can be optionally set by input argument

if nargin == 3
    c(6) = offset;
end

% data (nx2 double): first column is x positions and second is y positions
% x(i) and y(i) are coordinate pairs to evaluate Gaussian at
x = data(:,1);
y = data(:,2);

% out (column vector double): Gaussian function evaluated at x and y
out = c(1).*exp(-(x-c(2)).^2./(2*c(4)^2) - (y-c(3)).^2./(2*c(5)^2)) + c(6);

end