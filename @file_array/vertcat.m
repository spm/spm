function o = vertcat(varargin)
% Vertical concatenation of file_array objects.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: vertcat.m 7147 2017-08-03 14:07:01Z spm $


o = cat(1,varargin{:});
