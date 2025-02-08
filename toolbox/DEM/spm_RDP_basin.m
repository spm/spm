function [MDP,d,o,h,c] = spm_RDP_basin(MDP,S,chi)
% BMR of superordinate states based upon predicted outcomes
% FORMAT MDP = spm_RDP_basin(MDP,S,chi)
% MDP  - cell array of MDP structures
% S    - list of modality streams; e.g. [2,3]
% chi  - log prior for each modality; e.g., [-8,8]
%
% MDP  - reduced MDP
% d    - vector indicating whether deep states are transient states  [true]
% o    - vector indicating whether deep states are not orphan states [true]
% h    - list of goal states
% c    - list of cost states
%
% This routine implements a Bayesian model reduction in which the latent
% (generalised) states at the highest level of the model are retained only
% if they lie in the basins of attraction of generalised (deep) states
% specified by  intended states (i.e., indexed by hid) precluding
% transitions through costly states (indexed by cid).
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

% basins of atraction (i.e., paths to intended states under constraints)
%==========================================================================

% allowable transitions
%--------------------------------------------------------------------------
B      = sum(MDP{end}.b{1},3) > 0;
B(c,:) = false;

% contrained transitions under time reversal
%--------------------------------------------------------------------------
Nt     = 128;
Ns     = size(B,1);
P      = false(Nt,Ns);
P(1,h) = true;

for t = 1:Nt
    P(t + 1,:) =  any(B(P(t,:),:),1)';
end

% find states in basins of attraction
%--------------------------------------------------------------------------
R    = true(1,Ns);
R(c) = false;                            % not cost states         
R    = R & any(P);                       % and any path to hid

% restriction matrix for reduction of MDP
%--------------------------------------------------------------------------
j   = R;
R   = speye(Ns,Ns);
R   = R(:,j);
MDP = spm_RDP_compress(MDP,R);

% update contraints
%--------------------------------------------------------------------------
MDP = spm_set_goals(MDP,S,chi);
d   = any(sum(MDP{end}.b{1},3),1);
o   = any(sum(MDP{end}.b{1},3),2);
h   = MDP{end}.id.hid;
c   = MDP{end}.id.cid;
return
