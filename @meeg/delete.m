function this = delete(this)
% Delete files of an M/EEG dataset from disk and return unlinked object
% FORMAT this = delete(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if islinked(this)
    spm_unlink(fnamedat(this));
end
this = unlink(this);
spm_unlink(fullfile(this));
