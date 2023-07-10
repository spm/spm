function [I, P, t] = spm_mesh_ray_intersect(M, R)
% Compute the intersection of ray(s) and triangle(s)
% FORMAT [I, P, t] = spm_mesh_ray_intersect(M, R)
% M   - a GIfTI object or patch structure or numeric array [v1;v2;v3]
% R   - ray defined as a structure with fields 'orig' for origin and 'vec'
%       for direction, stored as column vectors
%
% I   - logical vector indicating intersection hit
% P   - coordinates of intersections [Mx3]
% t   - distance to hit triangles
%__________________________________________________________________________
%
% This function implements the Moller-Trumbore ray-triangle intersection
% algorithm:
% "Fast, Minimum Storage Ray-Triangle Intersection". Tomas Moller and Ben
% Trumbore (1997). Journal of Graphics Tools. 2: 21-28.
% https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
%__________________________________________________________________________
%
% M = gifti(fullfile(spm('Dir'),'canonical','scalp_2562.surf.gii'));
% R = struct('orig',[-100 100 -50]','vec',[150 -250 130]');
% [I,P] = spm_mesh_ray_intersect(M,R);
% spm_mesh_render(M);
% hold on
% p = plot3([R.orig(1) R.orig(1)+R.vec(1)],...
%     [R.orig(2) R.orig(2)+R.vec(2)],...
%     [R.orig(3) R.orig(3)+R.vec(3)],'-r','LineWidth',4);
% plot3(P(:,1),P(:,2),P(:,3),'*g','LineWidth',4);
%__________________________________________________________________________

% Guillaume Flandin
% Copyright (C) 2020-2022 Wellcome Centre for Human Neuroimaging


%-Default outputs
%--------------------------------------------------------------------------
P = zeros(0,3);
[t,u,v] = deal([]);

%-Options
%--------------------------------------------------------------------------
prec = 1e-7;

%-Face culling?
switch 'all'
    case 'all'
        culling = @(x) abs(x) >= prec;
    case 'front'
        culling = @(x) x >= prec;
    case 'back'
        culling = @(x) -x >= prec;
    otherwise
        error('Unknown culling option.');
end

%-Get vertices coordinates for all triangles
%--------------------------------------------------------------------------
if isa(M,'struct') || isa(M,'gifti')
    v1 = M.vertices(M.faces(:,1),:)';
    v2 = M.vertices(M.faces(:,2),:)';
    v3 = M.vertices(M.faces(:,3),:)';
else
    v1 = M(1:3,:);
    v2 = M(4:6,:);
    v3 = M(7:9,:);
end

%-Edges sharing vertex v1
%--------------------------------------------------------------------------
e1 = v2 - v1;
e2 = v3 - v1;

%-Computer determinant: if close to zero, ray lies in plane of triangle
%--------------------------------------------------------------------------
h = crossproduct(R.vec,e2);
d = dotproduct(e1,h);
I = culling(d);
if all(~I), return; end

%-Compute barycentric coordinate u and test bounds
%--------------------------------------------------------------------------
s = R.orig - v1;
u = dotproduct(s,h) ./ d;
I = I & (u >= 0) & (u <= 1);
if all(~I), return; end

%-Compute barycentric coordinate v and test bounds
%--------------------------------------------------------------------------
q = crossproduct(s,e1);
v = dotproduct(R.vec,q) ./ d;
I = I & (v >= 0) & (u + v <= 1);
if all(~I), return; end

%-Compute distance t and intersection point
%--------------------------------------------------------------------------
t = dotproduct(e2,q) ./ d;
I = I & (t >= prec);
if nargout > 1
    t = t(I);
    P = (R.orig + R.vec .* t)';
end


%==========================================================================
%-Cross product (using matrix multiplication)
%==========================================================================
function C = crossproduct(A,B)
if isequal(size(A),size(B))
    C = cross(A,B);
else
    if numel(A) == 3
        % skew-symmetric matrix
        a = [0 -A(3) A(2); A(3) 0 -A(1); -A(2) A(1) 0];
        C = a * B;
    elseif numel(B) == 3
        C = - crossproduct(B,A);
    else
        error('Number of vectors in A and B has to be 1-1, N-1, 1-N or N-N.');
    end
end


%==========================================================================
%-Dot product
%==========================================================================
function C = dotproduct(A,B)
if isequal(size(A),size(B))
    C = dot(A,B);
else
    if numel(A) == 3
        C = A' * B;
    elseif numel(B) == 3
        C = B' * A;
    else
        error('Number of vectors in A and B has to be 1-1, N-1, 1-N or N-N.');
    end
end
