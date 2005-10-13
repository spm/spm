function a = reshape(b,varargin)
% Overloaded reshape function for file_array objects
% _______________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

%
% $Id: reshape.m 253 2005-10-13 15:31:34Z guillaume $


if length(struct(b))~=1, error('Can only reshape simple file_array objects.'); end;

args = [];
for i=1:length(varargin),
	args = [args varargin{i}(:)'];
end;
if prod(args)~=prod(b.dim),
	error('To RESHAPE the number of elements must not change.');
end;
a = b;
a.dim = args;

