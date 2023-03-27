% Bayesian model reduction subroutine
%==========================================================================
function MDP = spm_MDP_motor_learning(MDP)
% FORMAT MDP = spm_MDP_motor_learning(MDP)
% MDP - generative model
%
% This routine illustrates a particular kind of structure learning, with a
% special focus on whether or not a particular path is controllable. This
% can be implemented simply and efficiently by evaluating the ELBO under
% different models of control (encoded by the indicator variables in
% MDP.k). And selecting the model with the greatest evidence. The a priori
% imperative for action is the expected information gain inherent in the
% expected free energy. If these prior beliefs about the consequences of
% action are realised, then the model evidence is high. If these prior
% beliefs are inconsistent with what the agent can realise, then model
% evidence will fall. The evidence in question here is the
% variiational free energy due to posterior beliefs over paths (MDP.Z).
% 
% In short, this routine updates structural knowledge about controllable
% dynamics by updating MDP.k.
% 
%
% See: spm_MDP_log_evidence.m, spm_MDP_VB_update and spm_MDP_VB_sleep.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_MDP_motor_learning.m 8439 2023-03-27 18:41:45Z guillaume $
%__________________________________________________________________________


% number of outsomes, states, controls and policies
%--------------------------------------------------------------------------
[Nf,Ns,Nu] = spm_MDP_size(MDP);
OPTIONS.A  = 0;                             % suppress explicit action
OPTIONS.B  = 0;                             % suppress backward pass
OPTIONS.N  = 0;                             % suppress neuronal responses
OPTIONS.G  = 1;                             % suppress graphics


% get combinations of controllable factors
%--------------------------------------------------------------------------
f      = find(Nu > 1);                      % potentially controllable factors
u      = spm_perm_mtx(numel(f));            % combinations of control
Nc     = size(u,1);                         % number of combinations
U      = zeros(Nc,Nf);                      % control indicator matrix
U(:,f) = u;

% motor babbling (evaluate the path integral of free energy  under
% different beliefs about controllable factors)
%==========================================================================
MDP.T  = 8;
MDP.s  = ones(Nf,1);
MDP.u  = ones(Nf,1);
for c  = 1:Nc
    mdp   = MDP;
    mdp.k = U(c,:);                         % what the agent thinks
    mdp   = spm_MDP_VB_XXX(mdp,OPTIONS);    % motor babbling
    F(c)  = sum(mdp.Z);                     % assess evidence: ELBO(u)

    % Graphics
    %----------------------------------------------------------------------
    if OPTIONS.G
        spm_figure('getwin','Motor learning'); clf
        spm_MDP_VB_trial(mdp);
    end
 
end

% select the best control model
%--------------------------------------------------------------------------
[F,c] = max(F);
MDP.k = U(c,:);

return






