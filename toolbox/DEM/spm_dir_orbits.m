function [u] = spm_dir_orbits(b,hid,P,Nx)
% information geometry of a likelihood mapping
% FORMAT [u] = spm_dir_orbits(b,hid,P,Nx)
% b      - transition tensor
% hid    - states to highlight
% P      - (Ns x Nt) sequence of (probabilistic) states
% Nx     - rotate to show first Nx (generalised) states
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

% Plausible transitions
%--------------------------------------------------------------------------
B     = b > 1/16;

% first orbit
%--------------------------------------------------------------------------
if nargin < 4, Nx = Ns; end
O       = false(Ns,1);
O(1:Nx) = true;
R       = diag(O + 1/2);

% state space (defined by graph Laplacian)
%--------------------------------------------------------------------------
b     = spm_detrend(b);
b     = b + b' + eye(size(b));
u     = spm_svd(R*b*R);

% G   = b + b';
% G   = G - diag(diag(G));
% G   = G - diag(sum(G));
% G   = expm(G);
% u   = spm_svd(G);

% Plot latent states
%--------------------------------------------------------------------------
plot(u(:,1),u(:,2),'o')

% And flow based upon transition probabilities (B)
%--------------------------------------------------------------------------
Xlim  = get(gca,'XLim');
Ylim  = get(gca,'YLim');
Pos   = get(gca,'Position');
for s = 1:Ns
    r = find(B(1:Nx,s));
    for i = 1:numel(r)

        P1   = [u(s,1), u(s,2)];
        P2   = [u(r(i),1), u(r(i),2)];

        X(1) = Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P1(1)-Xlim(1));
        X(2) = Pos(1)+(Pos(3))/(Xlim(2)-Xlim(1))*(P2(1)-Xlim(1));
        Y(1) = Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P1(2)-Ylim(1));
        Y(2) = Pos(2)+(Pos(4))/(Ylim(2)-Ylim(1))*(P2(2)-Ylim(1));

        annotation('arrow', X, Y)
    end
end
title('Orbits')

% Highlight states if necessary
%==========================================================================
if nargin > 1
    hold on
    plot(u(hid,1),u(hid,2),'.r','MarkerSize',32)
    hold off
end

% Add dynamics if specified
%==========================================================================
if nargin < 3, P = []; end
if numel(P)
    a     = axis;
    for t = 1:size(P,2)
        [~,s] = max(P(:,t));
        color = min(max([(1 - P(s,t)),1,(1 - P(s,t))],0),1);
        hold off, plot(u(hid,1),u(hid,2),'.r','MarkerSize',32), hold on
        plot(u(s,1),u(s,2),'.','MarkerSize',48,'Color',color)
        axis(a)
        M(t)  = getframe(gca);
    end

    % svae frames
    %----------------------------------------------------------------------
    h = gca;
    set(h(1),'Userdata',[])
    set(h(1),'Userdata',{M,8})
    set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

end

return