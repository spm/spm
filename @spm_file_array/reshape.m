function a = reshape(b,varargin)
% Overloaded reshape function for spm_file_array objects
% _______________________________________________________________________
% %W% John Ashburner %E%

args = [];
for i=1:length(varargin),
	args = [args varargin{i}(:)'];
end;
if prod(args)~=prod(b.dim),
	error('To RESHAPE the number of elements must not change.');
end;
a = b;
a.dim = args;

