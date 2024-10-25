function MDP = spm_daisy_chain(R,S,MDP,MISS)
% Selctive structure learning of a hierarchical POMDP
% FORMAT MDP = spm_daisy_chain(R,S,MDP,MISS)
% R    - cell array of outcomes
% S    - cell array (attracting set) of output sequences
% MDP  - cell array of (renormalising) models
% MISS - list of constrained (i.e. restricted) outcomes at lowest level
%
% MDP  - new MDP
% O    - cell array of new paths in MDP
%
% This routine uses fast structure learning to append new trajectories to a
% recursive or renormalising generative model specified as a cell array of
% Markov decision processes. It augments past (generalised) states and
% transitions with new states and transitions given new outcomes (R), under
% the constraint that each event ends on the attracting set (S)
%
% See: spm_fast_structure_learning.m
%      spm_faster_structure_learning.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_structure_learning.m 8454 2023-11-04 17:09:11Z karl $
%__________________________________________________________________________


% options for model inversion (and evaluation)
%==========================================================================
if nargin < 4, MISS = []; end

% start of attracting paths
%--------------------------------------------------------------------------
dcat  = @(O) fix(spm_cat([O(:,1:end - 1); O(:,2:end)])*4)';
U     = [];
for m = 1:numel(S)
    U(m,:) = dcat(S{m}(:,1:2));
end
V     = dcat(R);


% for each path in attracting set
%--------------------------------------------------------------------------
Nm    = numel(MDP);
Ne    = 2^(Nm - 1);
for m = find(ismember(U,V,'rows')')
    T     = find(ismember(V,U(m,:),'rows')');
    T     = T(T > Ne);
    for i = 1:numel(T)

        % pre-path: new path (N)
        %------------------------------------------------------------------
        j     = (1:Ne) - Ne + T(i) - 1 ;
        N     = R(:,j);

        % most likely outcomes
        %------------------------------------------------------------------
        o     = zeros(size(N));
        for t = 1:size(o,2)
            for g = 1:size(o,1)
                [~,k]  = max(N{g,t});
                o(g,t) = k;
            end
        end

        % and append
        %------------------------------------------------------------------
        if ~any(ismember(o',MISS','rows'))
            MDP = spm_merge_structure_learning([N,S{m}],MDP);
        end
    end
end

return



