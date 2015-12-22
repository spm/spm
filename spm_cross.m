function [Y] = spm_cross(X,x)
% Multidimensional cross (outer) product
% FORMAT [Y] = spm_cross(X,x)
%
% X  - numeric array
% x  - vector
%
% Y  - outer product
%
% See also: spm_dot
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cross.m 6654 2015-12-22 12:55:36Z spm $


% inner product
%--------------------------------------------------------------------------
if isvector(X), Y = X(:)*x(:)'; return, end

d   = size(X);
ind = repmat(':,',1,numel(d));
ind = ind(1:end - 1);

d(end + 1) = 1;
Y          = zeros(d);
for i = 1:numel(x)
    eval(['Y(' ind ',' num2str(i) ') = X*x(i);']);
end
