% Bayesian model reduction subroutine
%==========================================================================
function MDP = spm_MDP_motor_learning(MDP)
% FORMAT MDP = spm_MDP_motor_learning(MDP)
% MDP - generative model
% MDP.U - true control
%
% returns
% MDP.k - inferred control
%
% This routine illustrates a particular kind of structure learning, with a
% special focus on whether or not a particular path is controllable. This
% can be implemented simply and efficiently by evaluating the ELBO under
% different models of control (encoded by the indicator variables in
% MDP.k), and selecting the model with the greatest evidence. In this
% implementation, the outcomes are generated under true control (i.e.,
% under planning as inference). Using these outcomes, the evidence for
% various models of control (i.e., the indicator variables in MDP.k) is
% assessed. The evidence in question here is the variiational free energy
% due to posterior beliefs over paths (MDP.Z).
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

% priors over initial states and control (first state and control)
%--------------------------------------------------------------------------
MDP.s = ones(Nf,1);
MDP.u = ones(Nf,1);
for f = 1:Nf
    MDP.D{f} = sparse(1,1,1,Ns(f),1);        % initial state
    MDP.E{f} = sparse(1,1,1,Nu(f),1);        % initial control
end

% generate outcomes under true control
%--------------------------------------------------------------------------
MDP.T = 16;
mdp   = MDP;
mdp.k = MDP.U;                              % true contol
mdp   = spm_MDP_VB_XXX(mdp,OPTIONS);        % generate observations
o     = mdp.o;

% motor babbling (evaluate the path integral of free energy  under
% different beliefs about controllable factors)
%==========================================================================
for c  = 1:Nc

    mdp   = MDP;
    mdp.k = U(c,:);                         % what the agent thinks
    mdp.o = o;                              % what the agent sees
    mdp   = spm_MDP_VB_XXX(mdp,OPTIONS);    % motor babbling
    F(c)  = sum(mdp.F);                     % assess evidence: ELBO(u)

    % Graphics
    %----------------------------------------------------------------------
    if OPTIONS.G
        spm_figure('getwin','Motor learning'); clf
        spm_MDP_VB_trial(mdp);
    end
 
end

% select the best control model
%--------------------------------------------------------------------------
[f,c] = max(F);
MDP.k = U(c,:);

if OPTIONS.G
    spm_figure('getwin','Motor learning'); clf
    mdp   = spm_MDP_VB_XXX(MDP,OPTIONS);
    spm_MDP_VB_trial(mdp);

    subplot(3,2,3), imagesc(1 - U); title('Models of control')
    xlabel('hidden factor'); ylabel('model'); axis square
    subplot(3,2,4); bar(F - min(F)); title('Log evidence')
    xlabel('Models of control'); ylabel('ELBO (nats)'); axis square

end

return






