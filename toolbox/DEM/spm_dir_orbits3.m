function [u] = spm_dir_orbits3(b,hid)
% information geometry of a likelihood mapping
% FORMAT [u] = spm_dir_orbits3(b,hid)
% b      - transition tensor
% hid    - states to highlight
%
% This routine plots latent states in the first three principal
% eigenvectors of (roughly) the graph Laplacian of transitions among
% states. The transitions among states are supplied in terms of a
% transition tensor (that is summed over the third dimension or perhaps).
% This is a visual representation of orbits and paths in terms of follows
% in latent state space.
%
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging

% Approximate eigenvectors (u) of graph Laplacian (of transitions)
%==========================================================================
b  = spm_dir_norm(sum(b,3));

% state space (defined by graph Laplacian)
%--------------------------------------------------------------------------
b  = spm_detrend(b);
b  = b + b' + eye(size(b));
u  = spm_svd(b);

% Plot latent states
%--------------------------------------------------------------------------
hold off
plot3(u(:,1),u(:,2),u(:,3),'or'), hold on
plot3(u(:,1),u(:,2),u(:,3),':k'), hold on

% Highlight states if necessary
%==========================================================================
if nargin > 1
    plot3(u(hid,1),u(hid,2),u(hid,3),'.r','MarkerSize',32)
end
title('Orbits')
axis square

return