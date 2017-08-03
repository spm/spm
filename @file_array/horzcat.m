function o = horzcat(varargin)
% Horizontal concatenation of file_array objects
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

%
% $Id: horzcat.m 7147 2017-08-03 14:07:01Z spm $


o = cat(2,varargin{:});
