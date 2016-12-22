function this = rmdata(this)
% Deletes the data file and unlinks the header
% FORMAT this = rmdata(this)
% _________________________________________________________________________
% Copyright (C) 2011-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: rmdata.m 6976 2016-12-22 11:04:45Z vladimir $

if islinked(this)
    try
        delete(fnamedat(this));
    end
end

this = unlink(this);
