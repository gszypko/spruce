function [xbin,ybin,imgbin] = binImgForFit(x,y,img,numpts)
% x (1xn double): array of x values corresponding to columns of img
% y (1xm double): array of y values corresponding to rows of img
% img (mxn double): img(m,n) corresponds to x(n) and y(m)

xbin = x;
ybin = y;
imgbin = img;

while length(xbin) > numpts || length(ybin) > numpts
    % Ensure that x has an even number of elements
    if mod(length(xbin),2) ~= 0 
        xbin(end) = [];
        imgbin(:,end) = [];
    end

    % Ensure that y has an even number of elements
    if mod(length(ybin),2) ~= 0 
        ybin(end) = [];
        imgbin(end,:) = [];
    end

    if length(xbin) > numpts
        xbin = binMatrix(xbin,2,1);
        imgbin = binMatrix(imgbin,2,1);
    end

    if length(ybin) > numpts
        ybin = binMatrix(ybin,2,1);
        imgbin = binMatrix(imgbin,1,2);
    end
end

end
