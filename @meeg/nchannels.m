function res = nchannels(obj)
% Method for getting the number of channels
% FORMAT res = nchannels(obj)
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: nchannels.m 1125 2008-01-30 12:12:18Z vladimir $

res = length(obj.channels);