function res = delete(this)
% Delete the files of M/EEG dataset from the disk
% FORMAT res = delete(this)
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: save.m 2394 2008-10-23 15:38:38Z vladimir $

res = 1;

try
    delete(fullfile(path(this), fnamedat(this)));
    delete(fullfile(path(this), fname(this)));
catch
    res = 0;
end
