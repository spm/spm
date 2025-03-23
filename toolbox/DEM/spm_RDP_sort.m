function [MDP,j] = spm_RDP_sort(MDP,~)
% BMR using eigen-decomposition of allowable transitions 
% FORMAT [MDP,j] = spm_RDP_sort(MDP,['sort'])
% MDP  - hierarchal RDP
% sort - optional flag to supress pruning
% j    - sorted indices (at highest level)
%
% This routine uses the eigen-decomposition of allowable transitions (i.e.,
% flow among states) to identify steady-state distributions (i.e., the
% principal eigenvector with a real eigenvalue of unity). These
% nonequilibrium steady-state (NESS) distributions correspond to the
% attracting set; i.e., a pullback attractor.
%
%  Generalised states sorted according to the nonequilibrium steady-state
%  probability, which corresponds to the probability current through these
%  (generalised) states – because a path can never follow itself (and
%  therefore each generalised state will be vacated with probability one at
%  each time step).
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% eigenvalue decomposition of flows
%--------------------------------------------------------------------------
B     = spm_dir_norm(sum(MDP{end}.b{1},3) > 0);
[e,v] = eig(B,'nobalance');
Ns    = size(B,1);

% Nonequilibrium steady-state probabilities
%--------------------------------------------------------------------------
[~,j] = max(real(diag(v)));
p     = spm_dir_norm(abs(e(:,j)))';
j     = true(1,Ns);

% prune small NESS states if requested
%==========================================================================
if nargin < 2
    [~,k] = sort(p,'ascend');
    for i = k
        d     = j;
        d(i)  = false;
        if all(any(B(d,d)),1)
            j = d;
        end
    end
end

% Sort on NESS probability
%--------------------------------------------------------------------------
j     = find(j);
[~,k] = sort(p(j),'descend');
j     = j(k);

% restriction matrix (R) and reduction of MDP
%--------------------------------------------------------------------------
R   = speye(Ns,Ns);
MDP = spm_RDP_compress(MDP,R(:,j),'first');

return

