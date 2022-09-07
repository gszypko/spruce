function [] = move_ax(ax,dx,dy)

for i = 1:size(ax,1)
    for j = 1:size(ax,2)
        ax{i,j}.Position(1) = ax{i,j}.Position(1) + dx;
        ax{i,j}.Position(2) = ax{i,j}.Position(2) + dy;
    end
end

end