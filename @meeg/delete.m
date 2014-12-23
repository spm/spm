function this = delete(this)
% Delete files of an M/EEG dataset from disk and return unlinked object
% FORMAT this = delete(this)
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: delete.m 6293 2014-12-23 18:15:57Z guillaume $


if islinked(this)
    spm_unlink(fnamedat(this));
end
this = unlink(this);
spm_unlink(fullfile(this));
