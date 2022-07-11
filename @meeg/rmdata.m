function this = rmdata(this)
% Deletes the data file and unlinks the header
% FORMAT this = rmdata(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if islinked(this)
    try
        delete(fnamedat(this));
    end
end

this = unlink(this);
