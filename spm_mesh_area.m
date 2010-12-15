function A = spm_mesh_area(M,PF)
% Compute the surface area of a triangle mesh
% M        - patch structure: vertices and faces must be mx3 and nx3 arrays
%            or 3xm array of edge distances
% PF       - logical indicating whether to return surface area as a whole
%            or per face [default: false]
%
% A        - surface area
%__________________________________________________________________________
%
% Computed using numerically stable version of Heron's formula:
% See http://en.wikipedia.org/wiki/Heron%27s_formula
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin
% $Id: spm_mesh_area.m 4139 2010-12-15 18:31:49Z guillaume $

if isnumeric(M)
    A = M;
else
    A = M.vertices(M.faces',:);
    A = permute(reshape(A',3,3,[]),[2 1 3]);
    A = squeeze(sqrt(sum((A([1 2 3],:,:) - A([2 3 1],:,:)).^2,2)));
end
A = sort(A,1,'descend');
A = 1/4 * sqrt(( A(1,:) + ( A(2,:) + A(3,:) ) ) .* ...
               ( A(3,:) - ( A(1,:) - A(2,:) ) ) .* ...
               ( A(3,:) + ( A(1,:) - A(2,:) ) ) .* ...
               ( A(1,:) + ( A(2,:) - A(3,:) ) ));

if nargin < 2 || ~PF
    A = sum(A);
end
