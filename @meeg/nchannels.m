function res = nchannels(this)
% returns number of channels
% FORMAT res = nchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008-2012 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nchannels.m 5025 2012-10-31 14:44:13Z vladimir $

if this.montage.Mind==0
    res = length(this.channels);
else
    res = length(this.montage.M(this.montage.Mind).channels);
end