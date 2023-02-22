function [df] = derivative1D(dx,dy,f,axis)

df = zeros(size(f));
for ix = 2:size(f,2)-1 % x iterator
    for iy = 2:size(f,1)-1 % y iterator
        if strcmp(axis,'x')
            dfdx_left = (f(iy,ix) - f(iy,ix-1))/(dx(ix-1)/2+dx(ix)/2);
            dfdx_right = (f(iy,ix+1) - f(iy,ix))/(dx(ix)/2+dx(ix+1)/2);
            f_left = f(iy,ix) - dfdx_left*dx(ix)/2;
            f_right = f(iy,ix) + dfdx_right*dx(ix)/2;
            df(iy,ix) = (f_right/2 - f_left/2)/dx(ix);
        elseif strcmp(axis,'y')
            dfdy_left = (f(iy,ix) - f(iy-1,ix))/(dy(iy-1)/2+dy(iy)/2);
            dfdy_right = (f(iy+1,ix) - f(iy,ix))/(dy(iy)/2+dy(iy+1)/2);
            f_left = f(iy,ix) - dfdy_left*dy(iy)/2;
            f_right = f(iy,ix) + dfdy_right*dy(iy)/2;
            df(iy,ix) = (f_right/2 - f_left/2)/dy(iy);
        else
            error('<axis> must be <x> or <y>');
        end
    end
end


end