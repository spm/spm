function a = spm_file_array(fname,dim,dtype,scale)
% Function for creating spm_file_array objects.
% FORMAT a = spm_file_array(fname,dim,dtype, scale)
% a      - spm_file_array object
% fname  - filename
% dim    - dimensions
% dtype  - datatype (see spm_type)
% scale  - struct with scaling factors: scale.dim denotes the dimensions,
%          in which scaling takes place. scale.values contains a vector of
%          scaling values.
% _______________________________________________________________________
% @(#)spm_file_array.m	1.1 John Ashburner 04/06/28

if nargin < 3,
	error('Incorrect usage');
end;

dim = dim(:)';
dtype = dtype(1);

if nargin == 3
	scale = [];
end

if ~isempty(scale)
	for i = 1:length(scale.dim)
		if size(scale.values, i) ~= dim(scale.dim(i))
			error('Incorrect usage');
		end
	end
end

tmp = struct('fname',fname,'dim',dim,'dtype',dtype, 'scale', scale);
a = class(tmp, 'spm_file_array');
