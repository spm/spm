function res = nchannels(this)
% returns number of channels
% FORMAT res = nchannels(this)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nchannels.m 4432 2011-08-15 12:43:44Z christophe $

if this.montage.Mind==0
    res = length(this.channels);
else
    res = length(this.montage.M(this.montage.Mind).channels);
end