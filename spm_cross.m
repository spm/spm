function [Y] = spm_cross(X,x)
% Mulitimensional cross (outer) preoduct
% [Y] = spm_dot(X,x,DIM)
%
% X  - numeric array
% x  - vector
%
% Y  - outer product
%
% See also: spm_dot
%__________________________________________________________________________
% Copyright (C) 2005-2013 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cross.m 6606 2015-11-22 18:57:07Z karl $


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


