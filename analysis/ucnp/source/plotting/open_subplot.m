function [fig,ax,an,row,col] = open_subplot(num,name,value)

% open figure
if nargin == 3, fig = figure(name,value);
else, fig = figure;
end

% set figure position
[panel] = readSubplotPanels();
ind = find(num==panel.num);
fig.Position = panel.pos(ind,:);
fig.Color = [1 1 1];

% create axes
row = panel.row(ind);
col = panel.col(ind);
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
an.FontSize = 11;

end