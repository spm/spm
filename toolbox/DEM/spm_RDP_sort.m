function [MDP,B] = spm_RDP_sort(MDP)
% BMR using eigen-decomposition of allowable transitions 
% FORMAT [MDP,B] = spm_RDP_sort(MDP)
% MDP - hierarchal RDP
% B   - reduced transition priors (at highest level)
%
% This routine uses the eigen-decomposition of allowable transitions (i.e.,
% flow among states) to identify steady-state distributions (i.e., the
% principal eigenvector with a real eigenvalue of unity). These
% nonequilibrium steady-state (NESS) distributions correspond to the
% attracting set; i.e., a pullback attractor.
%
%  Generalised states that are not within the pullback attractor are
%  eliminated and the remaining states sorted according to the
%  nonequilibrium steady-state probability, which corresponds to the
%  probability current through these (generalised) states – because a path
%  can never  follow itself (and therefore each generalised state will be
%  vacated with probability one at each time step).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________

% flows at highest level
%--------------------------------------------------------------------------
n     = numel(MDP);
Ns    = size(MDP{n}.b{1},1);


% ensure there are no childless states
%==========================================================================
for i = 1:0

    % find states with children
    %----------------------------------------------------------------------
    B    = sum(MDP{n}.b{1},3);
    j    = any(B,1);

    % break if all states have children (i.e., NESS)
    %----------------------------------------------------------------------
    if sum(j) == Ns, break, end

    % restriction matrix (R) and reduction of MDP
    %----------------------------------------------------------------------
    Ns  = size(B,1);
    R   = speye(Ns,Ns);
    R   = R(:,j);
    MDP = spm_RDP_compress(MDP,R);

end

% eigenvalue decomposition of flows
%--------------------------------------------------------------------------
B     = spm_dir_norm(sum(MDP{n}.b{1},3));
[e,v] = eig(B,'nobalance');
v     = diag(v);

% Nonequilibrium steady-state  probabilities
%--------------------------------------------------------------------------
p     = [];
for i = find(real(v)' > (1 - exp(-16)))
    p(end + 1,:) = spm_dir_norm(abs(e(:,i)));
end

% Assign each state to one or more (n) attracting sets
%--------------------------------------------------------------------------
jp    = find(any(p > exp(-16),1));
[~,m] = max(p(:,jp),[],1);
jj    = [];
for n = 1:size(p,1)

    % Sort on NESS probability
    %----------------------------------------------------------------------
    jn    = jp(m == n);
    [~,k] = sort(p(n,jn),'descend');
    jn    = jn(k);
    pn    = p;
    for i = 1:numel(jn)

        % remove state if there is a reduced NESS solution
        %------------------------------------------------------------------
        [~,k] = min(pn(n,jn));
        j     = jn;
        j(k)  = [];
        if all(any(B(j,j),1))
            jn = j;
        else
            pn(n,jn(k)) = 1;
        end
    end
    jj    = [jj jn];
end

% restriction matrix (R) and reduction of MDP
%--------------------------------------------------------------------------
Ns  = size(B,1);
R   = speye(Ns,Ns);
R   = R(:,jj);
MDP = spm_RDP_compress(MDP,R);

return

