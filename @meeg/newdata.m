function obj = newdata(obj, fnamedat, dim, dtype)
% Initialises new data file_array
% FORMAT res = newdata(obj, fnamedat, [nchannels nsamples ntrials], dtype)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

% initialise new file_array object
d = file_array(fullfile(obj.path,...
    fnamedat),...
    dim,...
    dtype);

% physically initialise file
d(end,end,end) = 0;

% link into meeg object
obj.data.y = d;

% change filenames
obj.fname = [spm_str_manip(fnamedat, 'r') '.mat'];
obj.data.fnamedat = fnamedat;