function [u] = spm_dir_orbits(b,hid,N)
% information geometry of a likelihood mapping
% FORMAT [u] = spm_dir_orbits(b,hid,N)
% b      - transition tensor
% hid    - states to highlight
% P      - (Ns x Nt) sequence of (probabilistic) states
% N      - rotate to show first N (generalised) states
%
% This routine plots latent states in the first two principal eigenvectors
% of (roughly) the graph Laplacian of transitions among states. The
% transitions among states are supplied in terms of a transition tensor
% (that is summed over the third dimension or perhaps). The states are then
% connected with arrows denoting transitions over time. In effect, this is
% a visual representation of orbits and paths in terms of follows in latent
% state space.
%
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2012-2022 Wellcome Centre for Human Neuroimaging


% Approximate eigenvectors (u) of graph Laplacian (of transitions)
%==========================================================================
Ns    = size(b,1);
b     = spm_dir_norm(sum(b,3));
if nargin < 3, N = Ns; end

% Plausible transitions
%--------------------------------------------------------------------------
B     = b > 1/16;

% numerical approximation to NESS
%--------------------------------------------------------------------------
e     = mean(b^128,2);
p     = spm_dir_norm(e);

% eigenvalue decomposition of flows
%--------------------------------------------------------------------------
% [e,v] = eig(b,'nobalance');
% v     = diag(v);

% Nonequilibrium steady-state  probabilities
%--------------------------------------------------------------------------
% i = find(real(v)' > (1 - exp(-16)),1);
% p = spm_dir_norm(abs(e(:,i)));

[~,j] = sort(p,'descend');
j     = j(1:min(end - 2,N));
N     = numel(j);
B     = spm_dir_norm(B(j,j));

% remove hid states if necessary
%==========================================================================
if nargin > 1
    h      = false(1,Ns);
    h(hid) = true;
    h      = h(j);
    hid    = find(h);
end

% state space (defined by graph Laplacian)
%--------------------------------------------------------------------------
b     = B + B' + eye(size(B));
u     = spm_svd(b);

% Plot latent states
%--------------------------------------------------------------------------
plot3(u(:,1),u(:,2),u(:,3),'.r','MarkerSize',32), hold on

% And flow based upon transition probabilities (B)
%--------------------------------------------------------------------------
for j = 1:N
    r = find(B(:,j));
    for i = 1:numel(r)
        X   = [u(j,1), u(r(i),1)];
        Y   = [u(j,2), u(r(i),2)];
        Z   = [u(j,3), u(r(i),3)];
        line(X,Y,Z,'color','k')
    end
end
title('Orbits'), xlabel('1st dimension'), ylabel('2nd dimension'), axis square

% Highlight states if necessary
%==========================================================================
if nargin > 1
    hold on
    plot3(u(hid,1),u(hid,2),u(hid,3),'or','MarkerSize',16)
    hold off
end

return