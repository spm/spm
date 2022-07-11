function res = nchannels(this)
% returns number of channels
% FORMAT res = nchannels(this)
%__________________________________________________________________________

% Vladimir Litvak
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging


if this.montage.Mind==0
    res = length(this.channels);
else
    res = length(this.montage.M(this.montage.Mind).channels);
end
