function [img_filt,ind_x,ind_y] = filterHotPixels(img,kx,ky,px,py,threshold)

img_filt = img;
[img_sg,ind_x,ind_y] = sgfilt2D(img,kx,ky,px,py,false);

num_hot = 0;
for i = 1:numel(img)
    is_hot = abs(img(i) - img_sg(i)) > threshold*abs(img_sg(i));
    if is_hot
        num_hot = num_hot + 1;
        img_filt(i) = img_sg(i); 
    end
end
hot_percentage = num_hot/numel(img)*100;
disp(['Hot Pixels Removed: ' num2str(hot_percentage), '%'])

end