function [Mb] = binMatrix(M,dx,dy)
% M (mxn double): original matrix to be binned
% dx (1x1 double): number of elements to bin in x-direction
% dy (1x1 double): number of elements to bin in y-direction
% Mb (m/dx x m/dy): binned matrix

%% Function Notes
% In order to bin a vector (i.e., if M is a 1xn double), ensure that the number of elements to bin
% in the unused dimension (in this case, y) is one (i.e., dy = 1).

%% Bin the matrix M
[m,n] = size(M); % get sizes of original matrix

Mb = mean(reshape(M,dy,[]),1); % average y bins
Mb = reshape(Mb,m/dy,[]); % convert back into matrix

% Note that the transposes here are necessary
Mb = mean(reshape(Mb',dx,[]),1); % average x bins
Mb = reshape(Mb,n/dx,[])'; % convert back into matrix

end