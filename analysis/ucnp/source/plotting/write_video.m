function [] = write_video(filepath,frames)
% filepath (string): full path with filename but no extension
% frames (cell): an array of frames from getframe

vID = VideoWriter(filepath,'MPEG-4');

min_duration = 10;
min_framerate = length(frames)/min_duration;
if min_framerate < vID.FrameRate
    vID.FrameRate = min_framerate;
end

open(vID)
for i = 1:length(frames)
    writeVideo(vID,frames{i});
end
close(vID)

end