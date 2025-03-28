function [MDP,d,o,h,c,r] = spm_RDP_basin(MDP,S,chi,L)
% BMR of superordinate states based upon predicted outcomes
% FORMAT [MDP,d,o,h,c,r] = spm_RDP_basin(MDP,S,chi,L)
% MDP  - cell array of MDP structures
% S    - list of modality streams;    e.g., [2,3]
% chi  - log prior for each modality; e.g., [8,-8]
% L    - length of paths to and from goal states; e.g., [32,1] 
%
% MDP  - reduced MDP
% d    - vector indicating whether deep states have children [true]
% o    - vector indicating whether deep states have parents  [true]
% h    - list of goal states
% c    - list of cost states
% r    - radius of goal states
%
% This routine implements a Bayesian model reduction in which the latent
% (generalised) states at the highest level of the model are retained only
% if they lie in the basins of attraction of generalised (deep) states
% specified by intended states (i.e., indexed by hid) precluding
% transitions through costly states (i.e., indexed by cid).
%
% see: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
%__________________________________________________________________________
% Copyright (C) Karl Friston

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________

% get constraints
%--------------------------------------------------------------------------
MDP = spm_set_goals(MDP,S,chi);
h   = MDP{end}.id.hid;
c   = MDP{end}.id.cid;

% defaults
%--------------------------------------------------------------------------
if nargin < 4
    L = [32,1];
end

% basins of atraction (i.e., paths to intended states under constraints)
%==========================================================================

% preclude transitions to costly states
%--------------------------------------------------------------------------
B      = sum(MDP{end}.b{1},3) > 0;
B(c,:) = false;

% Parents of goal paths
%--------------------------------------------------------------------------
Nt     = max(L(1),1);
Ns     = size(B,1);
P      = false(Nt,Ns);
P(1,h) = true;
for t  = 1:Nt
    p  = any(B(P(t,:),:),1)';
    P(t + 1,:) = p;

    if ~any(p), break, end
end

% Children of goal paths
%--------------------------------------------------------------------------
Nt     = max(L(2),1);
C      = false(Nt,Ns);
C(1,h) = true;
for t  = 1:Nt
    p  = any(B(:,C(t,:)),2);
    C(t + 1,:) = p;
    
    if ~any(p), break, end
end

% remove costly paths
%--------------------------------------------------------------------------
R    = true(1,Ns);
R(c) = false;

% retain paths to or from goals with children
%--------------------------------------------------------------------------
R    = R & (any(P,1) | any(C,1));

% restriction matrix for reduction of MDP
%--------------------------------------------------------------------------
j   = R;
R   = speye(Ns,Ns);
MDP = spm_RDP_compress(MDP,R(:,j),'first');

% update contraints
%--------------------------------------------------------------------------
MDP = spm_set_goals(MDP,S,chi);

d   = any(sum(MDP{end}.b{1},3),1);       % with children
o   = any(sum(MDP{end}.b{1},3),2);       % with parents
h   = MDP{end}.id.hid;
c   = MDP{end}.id.cid;

% requested, get radii of goal paths
%--------------------------------------------------------------------------
if nargout > 5
    r = spm_RDP_radius(MDP,h,c);
end

return
