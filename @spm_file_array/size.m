function d = size(a,varargin)
% overloaded size function for spm_file_array objects.
% _______________________________________________________________________
% %W% John Ashburner %E%

if nargin==1,
	d = a.dim;
else,
	if varargin{1}<=length(a.dim),
		d = a.dim(varargin{1});
	else,
		d = 1;
	end;
end;
