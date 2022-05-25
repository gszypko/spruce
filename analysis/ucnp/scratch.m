clear
f = filesep;
path = 'E:\Grant-Gorman\data-mhd\05.23.22\set_0\grid-evol';
files = dir(path);
files([files.isdir]) = [];
num = extractAfter({files.name},'tEvol');
num = extractBefore(num,'.png');
num = str2double(num);
[~,sortind] = sort(num);
files = files(sortind);

v = VideoWriter('E:\Grant-Gorman\data-mhd\05.23.22\set_0\test.mp4','MPEG-4');
duration = 3;
v.FrameRate = (length(files))/duration;
open(v)
for j = 1:10
for i = 1:length(files)
    A = imread([path f files(i).name]);
    writeVideo(v,A)
end
end
close(v)
