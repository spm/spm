function res = delete(this)
% Delete the files of M/EEG dataset from the disk
% FORMAT res = delete(this)
%_______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: delete.m 5025 2012-10-31 14:44:13Z vladimir $

res = 1;

try
    if islinked(this)
        delete(fnamedat(this));
    end
    delete(fullfile(this));
catch
    res = 0;
end
