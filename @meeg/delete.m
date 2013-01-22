function this = delete(this)
% Delete the files of M/EEG dataset from the disk
% FORMAT this = delete(this)
% returns unlinked object
%_______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: delete.m 5200 2013-01-22 10:08:01Z vladimir $


if islinked(this)
    spm_unlink(fnamedat(this));
end
this = unlink(this);
spm_unlink(fullfile(this));

