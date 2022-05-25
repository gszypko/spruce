function [cax] = get_axis(fig,ax)

set(fig,'CurrentAxes',ax)
cax = gca;

end