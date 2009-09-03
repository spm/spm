function res = size(this, varargin)
% returns the dimensions of the data matrix
% FORMAT res = size(this, dim))
% _______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Vladimir Litvak
% $Id: size.m 3350 2009-09-03 13:19:20Z vladimir $

res = size(this.data.y, varargin{:});

if ntrials(this) == 1
    res = [res 1];
end
