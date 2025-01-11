function [MDP,B] = spm_RDP_sort(MDP,N)
% BMR using eigen-decomposition of allowable transitions 
% FORMAT [MDP,B] = spm_RDP_sort(MDP,N)
% MDP - hierarchal RDP
% N   - upper bound on number of states [default: Inf]
% B   - reduced transition priors (at highest level)
%
% This routine uses the eigen-decomposition of allowable transitions (i.e.,
% flow among states) to identify steady-state distributions (i.e., the
% principal eigenvector with a real eigenvalue of unity). These
% steady-state distributions correspond to the attracting set.
%
% If called without an upper bound on the number of latent states, this
% routine will simply reorder high-level states to foreground successive
% states that are most frequently visited at steady state.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% get predictive mapping from states to outcomes in trailing streams
%==========================================================================
if nargin < 2, N = Inf; end

% flows at highest level
%--------------------------------------------------------------------------
n     = numel(MDP);
B     = spm_dir_norm(sum(MDP{n}.b{1},3));

% eigenvalue decomposition of flows
%--------------------------------------------------------------------------
% [e,v] = eig(B,'nobalance');
% v     = diag(v);
% [~,j] = sort(real(v),'descend');
% p     = e(:,j).^2;
% [~,j] = sort(abs(p(:,1)),'descend');
% j     = j(1:min(end,N));



% numerical approximation to NESS
%--------------------------------------------------------------------------
e     = mean(B^128,2);
p     = spm_dir_norm(e.^2);
[~,j] = sort(p,'descend');
j     = j(1:min(end,N));

% re-order likelihood mappings from first stream
%--------------------------------------------------------------------------
Ns    = numel(j);
for g = 1:numel(MDP{n}.a)
    if MDP{n}.id.A{g} == 1
        MDP{n}.a{g} = MDP{n}.a{g}(:,j);
    end
end

% and transitition priors
%--------------------------------------------------------------------------
Nu    = size(MDP{n}.b{1},3);
B     = zeros(Ns,Ns,Nu);
for u = 1:Nu
    B(:,:,u) = MDP{n}.b{1}(j,j,u);
end
MDP{n}.b{1} = B;

return


