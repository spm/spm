function this = delete(this)
% Delete the files of M/EEG dataset from the disk
% FORMAT this = delete(this)
% returns unlinked object
%_______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: delete.m 5068 2012-11-19 15:00:07Z vladimir $


if islinked(this)
    delete(fnamedat(this));
end
this = unlink(this);
delete(fullfile(this));

