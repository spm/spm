function [R,g,f,r] = spm_RDP_radius(MDP,h,c)
% BMR of superordinate states based upon predicted outcomes
% FORMAT [R,g,f,r] = spm_RDP_radius(MDP,h,c)
% MDP  - cell array of MDP structures
% h    - list of goal states
% c    - list of cost states
%
% R    - Ajaceny matrix among h
% g    - goal states with children
% f    - goal states on attractor
% r    - radius
%
% This routine returns the radius of each intended (goal) state for the
% leading factor of a renormalising generative model. The radius is defined
% as the maximum number of transitions between an orphan state (with no
% parents) and the intended state. 
%
% In addition, it computes the adjacency among goal states in terms of
% whether one goal state can be accessed from another. This allows one to
% identify a subset of goal states that constitute an orbit (c.f., pullback
% attractor).
%
% see: https://en.wikipedia.org/wiki/Absorbing_Markov_chain
%__________________________________________________________________________
% Copyright (C) Karl Friston

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________

% basins of atraction (i.e., paths to intended states under constraints)
%==========================================================================

% preclude transitions to costly paths
%--------------------------------------------------------------------------
B      = sum(MDP{end}.b{1},3) > 0;
B(c,:) = false;

Nt    = 32;
Ns    = size(B,1);
Nh    = numel(h);
r     = zeros( 1,Nh);
g     = zeros( 1,Nh);
R     = zeros(Nh,Nh);
for i = 1:Nh

    % contrained transitions
    %--------------------------------------------------------------------------
    C         = false(Nt,Ns);
    C(1,h(i)) = true;

    for t  = 1:Nt
        p          = any(B(:,C(t,:)),2);
        C(t + 1,:) = p;

        if ~any(p), break, end
    end

    % contrained transitions under time reversal
    %----------------------------------------------------------------------
    P         = false(Nt,Ns);
    P(1,h(i)) = true;

    for t = 1:Nt
        p          = any(B(P(t,:),:),1)';
        P(t + 1,:) = p;

        if ~any(p), break, end
    end

    % radius (length of longest path to i-th goal)
    %----------------------------------------------------------------------
    r(i)   = find(any(P,2),1,'last');

    % path from other goals
    %----------------------------------------------------------------------
    R(i,:) = any(P(2:end,h),1);

    % path from other goals
    %----------------------------------------------------------------------
    g(i)   = any(any(C(2:end,h),1));

end

if nargout < 3, return, end

% number of goal states constituting an orbit
%--------------------------------------------------------------------------
S  = R^32;
f  = any(S,1) & any(S,2)';

return
