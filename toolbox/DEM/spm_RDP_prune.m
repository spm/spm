function [MDP,j] = spm_RDP_prune(MDP,hid)
% BMR of superordinate states based upon predicted outcomes
% FORMAT [MDP,j] = spm_RDP_prune(MDP,[hid])
% MDP  - cell array of MDP structures
% hid  - list of goal states
%
% MDP  - reduced MDP
% j    - list of deep states with children
%
% This routine prunes the leading factor of a renormalising generative
% model by removing all childless states recursively, until all states are
% transient (or there are no states).
%
% see: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
%__________________________________________________________________________
% Copyright (C) Karl Friston

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% basins of atraction (i.e., paths to intended states under constraints)
%==========================================================================

% remove childless states recursively 
%--------------------------------------------------------------------------
B     = sum(MDP{end}.b{1},3) > 0;
Ns    = size(B,1);
j     = 1:Ns;
for i = 1:Ns
    d = any(B,1);
    B = B(d,d);
    j = j(d);

    % break of all states have successors
    %----------------------------------------------------------------------
    if all(d), break, end
end

% remove states with a long (> Nt) path to goal states
%==========================================================================
if nargin > 1
    Nt    = 32;
    h     = sparse(1,hid,1,1,Ns);
    h     = h(j);
    P     = zeros(Nt,numel(j));
    for t = 1:Nt
        P(t,:) = h;
        h      = (h*B) > 0;
    end

    % retain paths to hid
    %----------------------------------------------------------------------
    j     = j(any(P,1));
end


% reduce MDP
%==========================================================================

% restriction matrix for reduction of MDP
%--------------------------------------------------------------------------
R   = speye(Ns,Ns);
MDP = spm_RDP_compress(MDP,R(:,j),'first');

return
