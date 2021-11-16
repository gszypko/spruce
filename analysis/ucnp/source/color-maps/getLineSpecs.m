function [out] = getLineSpecs(m)
%%
% m is the number of lines to be plotted on the same graph

s.col = plasmaModified(m); % each row of s.col contains RGB triplet
s.style = cell(m,1);        % line style for each line color (row of s.col)

ind = 0;
for i = 1:m
    ind = ind+1;
    styleInd = {'-'};
    markerInd = {'o'};

    if (m == 4 || m == 5)
        styleInd = {'-','--'};
        markerInd = {'o','s'};
        s.style{i} = styleInd{ind};
        s.marker{i} = markerInd{ind};
    elseif (m > 5 && m < 10)
        styleInd = {'-','--','-.'};
        markerInd = {'o','s','d'};
        s.style{i} = styleInd{ind};
        s.marker{i} = markerInd{ind};
    elseif (m >= 10)
        styleInd = {'-','--','-.',':'};
        markerInd = {'o','s','d','v'};
        s.style{i} = styleInd{ind};
        s.marker{i} = markerInd{ind};
    else
        s.style{i} = styleInd{ind};
        s.marker{i} = markerInd{ind};
    end

    if ind >= length(styleInd)
        ind = 0;
    end
end

for i = 1:length(s.style)
    out(i).style = s.style{i};
    out(i).marker = s.marker{i};
    out(i).col = s.col(i,:);
end

end
