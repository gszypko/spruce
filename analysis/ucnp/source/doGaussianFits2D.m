function [data] = doGaussianFits2D(data,~)

% determine which variable corresponds to ion density
if strcmp(data.config.eq_set,'ideal_2F')
    density_name = 'i_n';
else
    density_name = 'n';
end

% fit 2D Gaussian to ion density distribution for each time point
for i = 1:length(data.grids.time)    
    disp(['2D Gaussian Fits: ' num2str(i) '/' num2str(length(data.grids.time))])
    
    x = data.grids.x_vec;
    y = data.grids.y_vec;
    img = data.grids.vars(i).(density_name);
    fit = fitImgWithGaussian(x,y,img,0);
    fields = {'imgfit','fit','amp','x0','y0','sigx','sigy','offset'};
    for j = 1:length(fields)
        data.grids.gauss_fits(i).(fields{j}) = fit.(fields{j});
    end
end

end