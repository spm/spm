function BMA = DEM_demo_MDP_fit
% Demo of active inference for trust games
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
% inference (with variational Bayes) to model foraging for information in a
% three arm maze.  This demo illustrates the inversion of single-subject
% and group data to make inferences about subject-specific parameters -
% such as their prior beliefs about precision and utility. We first
% generate some synthetic data for a single subject and illustrate the
% recovery of key parameters using variational Laplace. We then consider
% the inversion of multiple trials from a group of subjects to illustrate
% the use of empirical Bayes in making inferences at the between-subject
% level. Finally, we demonstrate the use of Bayesian cross-validation to
% retrieve out-of-sample estimates (and classification of new subjects).
%
% In this example, an agent starts at the centre of a three way maze that
% is baited with a reward in one of the two upper arms. However, the
% rewarded arm changes from trial to trial.  Crucially, the agent can
% identify where the reward (US) is located by accessing a cue (CS) in the
% lower arm. This tells the agent whether the reward is on the left or the
% right upper arm. This means the optimal policy would first involve
% maximising information gain or epistemic value by moving to the lower arm
% and then claiming the reward thus signified. Here, there are eight hidden
% states (four locations times right or left reward), four control states
% (that take the agent to the four locations) and seven outcomes (three
% locations times two cues plus the centre).  The central location has an
% ambiguous or uninformative outcome, and the upper arms are rewarded
% probabilistically.
%
% see also: spm_MPD_VB.m, spm_dcm_mdp.m and spm_nlsi_Newton.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_fit.m 7679 2019-10-24 15:54:07Z spm $
 
% set up and preliminaries: first generate synthetic (single subject) data
%==========================================================================
rng('default')
  
% outcome probabilities: A
%--------------------------------------------------------------------------
% We start by specifying the probabilistic mapping from hidden states
% to outcomes.
%--------------------------------------------------------------------------
a      = .95;
b      = 1 - a;
A      = [1 1 0 0 0 0 0 0;    % ambiguous starting position (centre)
          0 0 a b 0 0 0 0;    % left arm selected and rewarded
          0 0 b a 0 0 0 0;    % left arm selected and not rewarded
          0 0 0 0 b a 0 0;    % right arm selected and not rewarded
          0 0 0 0 a b 0 0;    % right arm selected and rewarded
          0 0 0 0 0 0 1 0;    % informative cue - reward on right
          0 0 0 0 0 0 0 1];   % informative cue - reward on left
 
% controlled transitions: B{u}
%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions of hidden states
% under each action or control state. Here, there are four actions taking the
% agent directly to each of the four locations.
%--------------------------------------------------------------------------
B{1}  = [1 0 0 1; 0 1 0 0;0 0 1 0;0 0 0 0];
B{2}  = [0 0 0 0; 1 1 0 1;0 0 1 0;0 0 0 0];
B{3}  = [0 0 0 0; 0 1 0 0;1 0 1 1;0 0 0 0];
B{4}  = [0 0 0 0; 0 1 0 0;0 0 1 0;1 0 0 1];
for i = 1:4
    B{i} = kron(B{i},eye(2));
end
 
% priors: (utility) C
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of log
% probabilities. Here, the agent prefers rewarding outcomes
%--------------------------------------------------------------------------
c  = 2;
C  = [0 c -c c -c 0 0]';
 
% now specify prior beliefs about initial state, in terms of counts
%--------------------------------------------------------------------------
d  = kron([1 0 0 0],[1 1])';
 
% allowable policies (of depth T).  These are just sequences of actions
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  3  4  4  4  4
      1  2  3  4  2  3  1  2  3  4];
 
 
% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.V = V;                    % allowable policies
mdp.A = A;                    % observation model
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred states
mdp.D = d;                    % prior over initial states
mdp.s = 1;                    % initial state

mdp.alpha = 2;                % precision of action selection
mdp.beta  = 1;                % inverse precision of policy selection
 
% true parameters
%--------------------------------------------------------------------------
n         = 128;              % number of trials
i         = rand(1,n) > 1/2;  % randomise hidden states over trials
P.beta    = log(2);
P.C       = log(2);

MDP       = mdp;
MDP.C     = mdp.C*exp(P.C);
MDP.beta  = mdp.beta*exp(P.beta);


[MDP(1:n)] = deal(MDP);
[MDP(i).s] = deal(2);


% Solve to generate data
%==========================================================================
MDP  = spm_MDP_VB(MDP);
 
% illustrate behavioural responses - single trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1a'); clf
spm_MDP_VB_trial(MDP(1));
 
% illustrate behavioural responses and neuronal correlates over trials
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1b'); clf
spm_MDP_VB_game(MDP);

%--------------------------------------------------------------------------
% This completes the generation of data. We now turn to the estimation of
% subject specific preferences and precision encoded by the parameters
% beta and C. Model parameters here are log scaling parameters that allow
% for increases or decreases in the default prior values.
%--------------------------------------------------------------------------


% Invert to recover parameters (preferences and precision)
%==========================================================================
DCM.MDP   = mdp;                  % MDP model
DCM.field = {'beta','C'};         % parameter (field) names to optimise
DCM.U     = {MDP.o};              % trial specification (stimuli)
DCM.Y     = {MDP.u};              % responses (action)

DCM       = spm_dcm_mdp(DCM);

% compare true values with posterior estimates
%--------------------------------------------------------------------------
subplot(2,2,4),hold on
bar(spm_vec(P),1/4)
set(gca,'XTickLabel',DCM.field)
set(gcf,'Name','Figure 2','Tag','Figure 2')


% now repeat using subsets of trials to illustrate effects on estimators
%==========================================================================
DCM.field = {'beta'};
n         = [8 16 32 64 128];
for i = 1:length(n)
    DCM.U = {MDP(1:n(i)).o};
    DCM.Y = {MDP(1:n(i)).u};
    DCM   = spm_dcm_mdp(DCM);
    Ep(i,1) = DCM.Ep.beta;
    Cp(i,1) = DCM.Cp;
end

% plus results
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
subplot(2,1,1), spm_plot_ci(Ep(:),Cp(:)), hold on
plot(1:length(n),(n - n) + P.beta,'k'),       hold off
set(gca,'XTickLabel',n)
xlabel('number of trials','FontSize',12)
ylabel('conditional estimate','FontSize',12)
title('Dependency on trial number','FontSize',16)
axis square


% now repeat but over multiple subjects with different betsa
%==========================================================================

% generate data and a between subject model with two groups of eight
% subjects
%--------------------------------------------------------------------------
N     = 8;                             % numbers of subjects per group
X     = kron([1 1;1 -1],ones(N,1));    % design matrix
h     = 4;                             % between subject log precision
n     = 128;                           % number of trials
i     = rand(1,n) > 1/2;               % randomise hidden states 

clear MDP
[MDP(1:n)] = deal(mdp);
[MDP(i).s] = deal(2);

for i = 1:size(X,1)
    
    % true parameters - with a group difference of one quarter
    %----------------------------------------------------------------------
    beta(i)    = X(i,:)*[0; 1/4] + exp(-h/2)*randn;
    [MDP.beta] = deal(exp(beta(i)));

    % solve to generate data
    %----------------------------------------------------------------------
    DDP        = spm_MDP_VB(MDP);      % realisation for this subject
    DCM.U      = {DDP.o};              % trial specification (stimuli)
    DCM.Y      = {DDP.u};              % responses (action)
    GCM{i,1}   = DCM;
    
    % plot behavioural responses
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 4'); clf
    spm_MDP_VB_game(DDP);drawnow
    
end


% Bayesian model inversion 
%==========================================================================
GCM  = spm_dcm_fit(GCM);

% plot subject specific estimates and true values
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4');
subplot(3,1,3)
for i = 1:length(GCM)
    qP(i) = GCM{i}.Ep.beta;
end
plot(beta,beta,':b',beta,qP,'.b','MarkerSize',32)
xlabel('true parameter','FontSize',12)
ylabel('conditional estimate','FontSize',12)
title('Subject specific estimates','FontSize',16)
axis square

                    
% hierarchical (empirical) Bayes
%==========================================================================

% second level model
%--------------------------------------------------------------------------
M    = struct('X',X);

% BMA - (second level)
%--------------------------------------------------------------------------
PEB  = spm_dcm_peb(GCM,M);
BMA  = spm_dcm_peb_bmc(PEB);

subplot(3,2,4),hold on, bar(1,1/4,1/4), set(gca,'XTickLabel',DCM.field)
subplot(3,2,2),hold on, bar(1,1/4,1/4), set(gca,'XTickLabel',DCM.field)


% posterior predictive density and cross validation
%==========================================================================
spm_dcm_loo(GCM,M,DCM.field);








