function [x,dx] = getNonUniformGrids(grid_number,domain_size,growth,spread,is_symmetric)
% grid_number (int): length of x (i.e., number of grids)
% domain_size (double): sets the limits of the domain in whatever unit you desire
% growth (double): sets the growth rate of the grid size
% spread (double): controls how far from the origin the grids grow appreciably
% is_symmetric (boolean): if true, reflect grids across the origin to create a symmetric vector

% Description:
%   - This function creates corresponding grid positions and sizes for a one dimensional vector.
%   - The grid size always grows with distance from the origin.
%   - The grid size is always scaled so that the total domain size is correct.

% default arguments
if nargin < 5, is_symmetric = true; end
if nargin < 4, spread = 0.9999; end
if nargin < 3, growth = 1.05; end

% iterate through growth equation
iterations = grid_number-1;
if is_symmetric, iterations = ceil(grid_number/2); end
x = zeros(1,iterations);
dx = zeros(1,iterations);
dx(1) = 1; % default start value for the growth equation, to be scaled to appropriate units later
for i = 2:iterations
    dx(i) = 1 + growth*(dx(i-1)-spread);
    x(i) = x(i-1)+dx(i-1)/2+dx(i)/2;
end

% scale to appropriate units
scale_fac = domain_size./max(x);
dx = dx.*scale_fac;
x = x.*scale_fac;

% account for whether domain is symmetric or not
if is_symmetric
    dx = [dx(end:-1:2) dx];
    x = [-x(end:-1:2) x];
end

end