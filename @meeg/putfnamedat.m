function obj = putfnamedat(obj, fnamedat)
% Method for putting a new file name for the *.dat file
% FORMAT obj = putfnamedat(obj, fnamedat)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

obj.data.y.fname = fullfile(obj.path, fnamedat);
