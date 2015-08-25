function M = spm_mesh_inflate(M,T,S)
% Surface mesh inflation
% FORMAT M = spm_mesh_inflate(M,T,S)
%
% M        - surface mesh structure (see patch) or GIfTI object
%            or handle to a patch in a figure
% T        - number of time steps [default: Inf (auto)]
% S        - update display every S time steps [default: 0 (never)]
%__________________________________________________________________________
% Copyright (C) 2009-2011 Wellcome Trust Centre for Neuroimaging

% Guillaume Flandin & Jean Daunizeau
% $Id: spm_mesh_inflate.m 6535 2015-08-25 11:45:26Z vladimir $


if nargin < 3, S = 0; end

if ishandle(M)
    v = double(get(M,'Vertices'));
    f = get(M,'Faces');
    p = M;
else
    v = double(M.vertices);
    f = double(M.faces);
    
    h = figure('visible','off');
    a = axes('parent',h);
    if isa(M,'gifti')
        p = patch(export(M,'patch'),'parent',a,'visible','off');
    else
        p = patch(M,'parent',a,'visible','off');
    end
    S = 0;
end

% Parameters
%--------------------------------------------------------------------------
b = 0.5;
w = 0.05;

if nargin < 2 || isinf(T)
    T = floor(size(v,1) * 0.003 - 2);
    T = max(T,1);
end

% Compute (normalised) adjacency matrix
%--------------------------------------------------------------------------
A = spm_mesh_adjacency(f);
A = sparse(1:size(v,1),1:size(v,1),1./sum(A,2)) * A;

% Compute bounding box
%--------------------------------------------------------------------------
minxyz = min(v);
maxxyz = max(v);

% Iteratively apply forces to vertices
%--------------------------------------------------------------------------
for i=1:T
    
    % Compute unit normals
    %----------------------------------------------------------------------
    N = spm_mesh_normals(p,1);

    % Compute smoothing force
    %----------------------------------------------------------------------
    mv = A*v - v;

    % Update vertices position
    %----------------------------------------------------------------------
    v = v + b * (w*repmat(sum(mv.*N,2),1,3).*N + (1-w)*mv);
    v = mean((maxxyz - minxyz)./(max(v) - min(v))) * v;
    
    set(p,'Vertices',v);
    
    % Update display
    %----------------------------------------------------------------------
    if ~mod(i,S)
        axis(get(p,'parent'),'image');
        drawnow
    end
    
end

% Cleanup
%--------------------------------------------------------------------------
if ~ishandle(M)
    M.vertices = v;
    close(h);
end
