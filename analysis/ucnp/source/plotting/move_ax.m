function [] = move_ax(ax,dx,dy)

for i = 1:length(ax)
    ax{i}.Position(1) = ax{i}.Position(1) - dx;
    ax{i}.Position(2) = ax{i}.Position(2) - dy;
end

end