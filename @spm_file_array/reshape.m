function a = reshape(b,varargin)
% Overloaded reshape function for spm_file_array objects
% _______________________________________________________________________
% @(#)reshape.m	1.1 John Ashburner 04/06/28

args = [];
for i=1:length(varargin),
	args = [args varargin{i}(:)'];
end;
if prod(args)~=prod(b.dim),
	error('To RESHAPE the number of elements must not change.');
end;
a = b;
a.dim = args;

