function this = rmdata(this)
% Deletes the data file and unlinks the header
% FORMAT this = rmdata(this)
% _________________________________________________________________________
% Copyright (C) 2011-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: rmdata.m 5025 2012-10-31 14:44:13Z vladimir $

if islinked(this)
    try
        delete(fullfile(path(this), fnamedat(this)));
    end
end

this = unlink(this);
