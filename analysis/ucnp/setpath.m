function [gitdir] = setpath(opt)
    % update the Matlab search path to include the required subfolders of the LIF analysis code 
    % repository
    
    if nargin < 1, opt = true; end
    
    f = filesep;
    fullpath = mfilename('fullpath');
    gitdir = extractBefore(fullpath,[f mfilename]);
    folders = {'','file-io','img-fits','other','plasma-quantities','plotting','ion-holes'};
    if opt == true
        addpath(gitdir)
        for i = 1:length(folders)
            addpath([gitdir f 'source' f folders{i}])
        end
    end
    
    end