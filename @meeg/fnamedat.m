function res = fnamedat(obj)
% Method for getting file name of data file
% FORMAT res = fnamedat(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel

res = spm_str_manip(obj.data.y.fname, 't');
