function [u] = spm_merge(p,q);
% image merging {c.f transparency overlay}
% FORMAT [u] = spm_merge(p,q);
% p - matrix {m x n}
% q - matrix {m x n}
% u - output {2m x 2n}
%___________________________________________________________________________
%
% spm_merge merges two images [matrices] in working memory,  The first
% image is considered as a 'foreground' image and if the voxel values
% are not zero then these voxel values are interspersed with those of
% the second, background image.
%
% The resulting image is usually displayed with the split {128,3} colormap
%
%__________________________________________________________________________
% %W% %E%

% interpolate images to twice their size
%---------------------------------------------------------------------------
[m n]      = size(p);
i          = max([0 min(p(p > 0))]);
j          = max([0 min(q(q > 0))]);
p          = spm_resize(p,2*m,2*n);
q          = spm_resize(q,2*m,2*n);
p          = p.*(p > i);
q          = q.*(q > j);

% create indices
%---------------------------------------------------------------------------
x          = zeros(1,2*m);
y          = zeros(1,2*n);
x([1:m]*2) = ones(1,m);
y([1:n]*2) = ones(1,n);
d          = toeplitz(x,y);

% combine where q > 0
%---------------------------------------------------------------------------
d          = find(d & p);
q(d)       = p(d);
u          = reshape(q,2*m,2*n);
