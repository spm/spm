function d = size(a,varargin)
% overloaded size function for spm_file_array objects.
% _______________________________________________________________________
% @(#)size.m	1.1 John Ashburner 04/06/28

if nargin==1,
	d = a.dim;
else,
	if varargin{1}<=length(a.dim),
		d = a.dim(varargin{1});
	else,
		d = 1;
	end;
end;
