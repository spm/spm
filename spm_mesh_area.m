function A = spm_mesh_area(M,P)
% Compute the surface area of a triangle mesh
% FORMAT A = spm_mesh_area(M,P)
% M        - patch structure: vertices and faces must be mx3 and nx3 arrays
%            or 3xm array of edge distances
% P        - return overall surface area, or per face, or per vertex
%            one of {'sum','face','vertex'} [default: 'sum']
%
% A        - surface area
%__________________________________________________________________________
%
% Computed using numerically stable version of Heron's formula:
% See https://www.wikipedia.org/wiki/Heron%27s_formula
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2010-2022 Wellcome Centre for Human Neuroimaging


if nargin < 2
    P = 'sum';
elseif islogical(P)
    % Backward compatibility with previous syntax
    if P
        P = 'face';
    else
        P = 'sum';
    end
end

if isnumeric(M)
    A = M;
else
    A = M.vertices(M.faces',:);
    A = permute(reshape(A',3,3,[]),[2 1 3]);
    A = squeeze(sqrt(sum((A([1 2 3],:,:) - A([2 3 1],:,:)).^2,2)));
end

A = sort(A,1,'descend');
A = ( A(1,:) + ( A(2,:) + A(3,:) ) ) .* ...
    ( A(3,:) - ( A(1,:) - A(2,:) ) ) .* ...
    ( A(3,:) + ( A(1,:) - A(2,:) ) ) .* ...
    ( A(1,:) + ( A(2,:) - A(3,:) ) );
A(A<0) = 0;
A = 1/4 * sqrt(A);

if strcmp(P,'sum')
    A = sum(A);
elseif strcmp(P,'face')
    A = A;
elseif strcmp(P,'vertex')
    AV = zeros(size(M.vertices,1),1);
    for i=1:3
        AV = AV + accumarray(M.faces(:,i),A,size(AV));
    end
    A = AV / 3; % each triangle contributes a third to each of its vertices
else
    error('Unknown option.');
end
