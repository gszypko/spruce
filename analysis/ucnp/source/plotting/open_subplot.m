function [fig,ax,an] = open_subplot(row,col,num)

if nargin < 3, num = row*col; end
if num > row*col, error('<num> must be smaller than <row>*<col>.'); end

% open figure
fig = figure;
fig.Position = [471.6667  251.6667  970.0000  672.0000];
fig.Color = [1 1 1];

% create axes
ax = cell(row,col);
iter = 0;
for i = 1:row
    for j = 1:col
        if iter > num - 1, break, end
        iter = iter + 1;
        ax{i,j} = subplot(row,col,iter);
    end
end

% create annotation
an = annotation('textbox');
an.Position = [0.1595    0.9514    0.7230    0.0371];
an.HorizontalAlignment = 'center';
an.VerticalAlignment = 'middle';
an.LineStyle = 'none';
an.FontSize = 12;

end