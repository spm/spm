function d = size(a,varargin)
% overloaded size function for spm_file_array objects.
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% John Ashburner
% $Id: size.m 253 2005-10-13 15:31:34Z guillaume $

if nargin==1,
	d = a.dim;
else,
	if varargin{1}<=length(a.dim),
		d = a.dim(varargin{1});
	else,
		d = 1;
	end;
end;
