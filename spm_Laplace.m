function u = spm_Laplace(u)
% Solve Laplace's equation on a regular grid
% FORMAT u = spm_Laplace(u)
% u        - potential field as a 3D array with values:
%            Inf: interior points (unknown values)
%            NaN: insulated boundaries
%            <b>: Dirichlet boundary conditions
%
% u        - filled-in potential field using Laplace's equation
%__________________________________________________________________________
%
% Potential field u should not have unknown values (Inf) at the first order
% boundary of the 3D array. Set them as insulated boundaries (NaN) if
% needed.
%
% See:
%
% Laplace's Equation in 2 and 3 Dimensions
% Douglas Wilhelm Harder, University of Waterloo, Canada
% https://ece.uwaterloo.ca/~ne217/Laboratories/05/5.LaplacesEquation.pptx
%__________________________________________________________________________

% Nicole Labra Avila, Guillaume Flandin
% Copyright (C) 2023 Wellcome Centre for Human Neuroimaging


%-Find interior points (labelled Inf)
%--------------------------------------------------------------------------
ip = find(isinf(u));
N  = numel(ip);

%-Build seven-point stencil layout [N (x 3) x 6]
%--------------------------------------------------------------------------
[x, y, z]    = ind2sub(size(u), ip);
XYZ          = repmat([x, y, z], [1, 1, 6]);
XYZ(:, 1, 1) = XYZ(:, 1, 1) - 1;
XYZ(:, 1, 2) = XYZ(:, 1, 2) + 1;
XYZ(:, 2, 3) = XYZ(:, 2, 3) - 1;
XYZ(:, 2, 4) = XYZ(:, 2, 4) + 1;
XYZ(:, 3, 5) = XYZ(:, 3, 5) - 1;
XYZ(:, 3, 6) = XYZ(:, 3, 6) + 1;
XYZ = squeeze(sub2ind(size(u), XYZ(:, 1, :), XYZ(:, 2, :), XYZ(:, 3, :)));

%-Types of neighbouring points [N x 6]
%--------------------------------------------------------------------------
t     = u(XYZ);

%-Indices of unknown neighbouring points
%--------------------------------------------------------------------------
w     = zeros(size(u));
w(ip) = 1:N;
m     = w(XYZ);
idx   = m > 0;
m     = m(idx);    % column indices
[k,~] = find(idx); % row indices

%-Build Laplacian matrix L [N x N]
%--------------------------------------------------------------------------
i = [(1:N)'; k];
j = [(1:N)'; m];
v = [-(6 - sum(isnan(t), 2)); ones(numel(k), 1)];
L = sparse(i, j, v, N, N);

%-Build right-hand side b (Dirichlet boundary condition) [N x 1]
%--------------------------------------------------------------------------
t(~isfinite(t)) = 0;
b = -sum(t, 2);  

%-Solve Laplace's equation Lu=b
%--------------------------------------------------------------------------
u(ip) = L \ b;
