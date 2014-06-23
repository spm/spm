function MDP = DEM_demo_MDP_maze
% Demo of active inference for trust games
%__________________________________________________________________________
%
% This routine uses the Markov decision process formulation of active
% inference (with variational Bayes) to model foraging for information in a
% three arm maze.  This demo illustrates variational free energy
% minimisation in the context of Markov decision processes, where the agent
% is equipped with prior beliefs that it will minimise expected free energy
% in the future. This free energy is the free energy of future sensory
% states expected under the posterior predictive distribution. It can be
% regarded as a generalisation of the variational formulation of KL control
% in which information gain or epistemic value is formulated explicitly.
%
% In this example, the agent starts at the centre of a three way maze
% which is baited with a reward in one of the two upper arms. However, the
% rewarded arm changes from trial to trial.  Crucially, the agent can
% identify where the reward (US) is located by accessing a cue (CS) in the
% lower arm. This tells the agent whether the reward is on the left or the
% right upper arm.  This means the optimal policy would first involve
% maximising information gain or epistemic value by moving to the lower arm
% and then claiming the reward this signified. Here, there are eight hidden
% states (four locations times right or left reward), four control states
% (that take the agent to the four locations) and 16 outcomes (four
% locations times two cues times two rewards).  The central location has an
% ambiguous or uninformative cue outcome, while the upper arms are rewarded
% probabilistically with an 80% schedule.
%
% A single trial is simulated followed by an examination of dopaminergic
% responses to conditioned and unconditioned stimuli (cues and rewards). A
% hierarchical version is then implemented, in which the mapping between
% locations in the generative model and the generative process is unknown
% and has to be learned.
%
% see also: spm_MPD_game
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_maze.m 6064 2014-06-23 09:39:46Z karl $

% set up and preliminaries
%==========================================================================
% rng('default')

% observation probabilities
%--------------------------------------------------------------------------
a      = .9;
A{1,1} = [.5 .5; .5 .5; 0 0; 0 0];
A{2,2} = [0 0; 0 0; a (1 - a); (1 - a) a];
A{3,3} = [0 0; 0 0; (1 - a) a; a (1 - a)];
A{4,4} = [1 0; 0 1; 0 0; 0 0];
A      = spm_cat(A);

% transition probabilities
%--------------------------------------------------------------------------
for i = 1:4
    B{i} = zeros(4,4);
    B{i}([2 3],[2 3]) = eye(2);
    B{i}(i,[1 4])     = 1;
    B{i} = kron(B{i},eye(2));
end

% priors: softmax(utility)
%--------------------------------------------------------------------------
c  = 2;
C  = spm_softmax(kron(ones(4,1),[0; 0; c; -c]));

% prior beliefs about initial state
%--------------------------------------------------------------------------
D  = kron([1 0 0 0],[1 1]/2)';

% true initial state
%--------------------------------------------------------------------------
S  = kron([1 0 0 0],[1 0])';


% allowable policies (of depth T)
%--------------------------------------------------------------------------
V  = [1  1  1  1  2  2  2  2  3  3  3  3  4  4  4  4
      1  2  3  4  1  2  3  4  1  2  3  4  1  2  3  4
      1  2  3  4  1  2  3  4  1  2  3  4  1  2  3  4];

 
% MDP Structure
%==========================================================================
MDP.N = 8;                          % number of variational iterations
MDP.S = S;                          % true initial state
MDP.A = A;                          % observation model
MDP.B = B;                          % transition probabilities (priors)
MDP.C = C;                          % terminal cost probabilities (priors)
MDP.D = D;                          % initial state probabilities (priors)
MDP.V = V;                          % allowable policies

MDP.alpha = 64;                     % gamma hyperparameters
MDP.beta  = 4;

% Solve - an example game
%==========================================================================
spm_figure('GetWin','Figure 1'); clf
MDP.plot = gcf;
MDP      = spm_MDP_game(MDP,'FE');


% return

% different formulations of optimality as a function of preference
%==========================================================================
spm_figure('GetWin','Figure 2'); clf
MDP.plot = 0;
MDP.N    = 4;

c     = linspace(0,2,8);
d     = kron(ones(4,1),[0; 0; 1; 0]);
for i = 1:length(c)
    
    % preference
    %----------------------------------------------------------------------
    MDP.C = spm_softmax(kron(ones(4,1),[0; 0; c(i); -c(i)]));
    
    
    % simulate trials and record outcomes
    %----------------------------------------------------------------------
    for j = 1:32
        
        s       = rand > 1/2;
        MDP.S   = kron([1 0 0 0],[s ~s])';
        
        MDP     = spm_MDP_game(MDP,'FE');
        FE(j,i) = d'*MDP.O(:,end);
        MDP     = spm_MDP_game(MDP,'KL');
        KL(j,i) = d'*MDP.O(:,end);
        MDP     = spm_MDP_game(MDP,'RL');
        RL(j,i) = d'*MDP.O(:,end);
    end
    
end
MDP.S  = S;
MDP.C  = C;

% posterior beliefs about hidden states (prosocial versus nonsocial)
%--------------------------------------------------------------------------
subplot(2,1,1)
bar(c,[mean(FE); mean(KL); mean(RL)]'*100), hold on
plot(c,c*0 + 100/3,'--k'), hold on
plot(c,c*0 + 100*a,'-.k'), hold off
title('Performance','FontSize',16)
xlabel('Prior preference','FontSize',12)
ylabel('success rate (%)','FontSize',12)
spm_axis tight, axis square
legend({'FE','KL','RL'})


% dopamine responses to US and CS
%==========================================================================
MDP.N = 16;
MDP.a = [4 2 2];
MDP.o = [1 ((4 - 1)*4 + 1) ((2 - 1)*4 + 3)];
MDP   = spm_MDP_game(MDP,'FE');

% axis
%--------------------------------------------------------------------------
ax    = [0.5*MDP.N 2.5*MDP.N 0 4];

spm_figure('GetWin','Figure 2');
subplot(2,1,2)
plot(MDP.da,'k'), hold on

MDP.a = [1 2 2];
MDP.o = [1 ((1 - 1)*4 + 1) ((2 - 1)*4 + 3)];
MDP   = spm_MDP_game(MDP,'FE');

spm_figure('GetWin','Figure 2');
subplot(2,1,2)
plot(MDP.da,'r'), hold off
title('Dopamine responses','FontSize',16)
xlabel('Peristimulus time','FontSize',12)
ylabel('Response','FontSize',12)
axis square, axis(ax)


% different levels of priors and uncertainty
%==========================================================================
spm_figure('GetWin','Figure 3'); clf; 
subplot(2,1,1)
MDP.a = [4 2 2];
MDP.o = [1 ((4 - 1)*4 + 1) ((2 - 1)*4 + 3)];
c     = linspace(0,2,4);
for i = 1:length(c)
    
    % preference
    %----------------------------------------------------------------------
    MDP.C = spm_softmax(kron(ones(4,1),[0; 0; c(i); -c(i)]));
    
    % simulate trials and record outcomes
    %----------------------------------------------------------------------
    MDP  = spm_MDP_game(MDP,'FE');
    col  = [1 1 1]*(length(c) - i)/length(c);
    plot(MDP.da,'Color',col), hold on

end
title('Preference (utility)','FontSize',16)
xlabel('Peristimulus time','FontSize',12)
ylabel('Response','FontSize',12)
axis square, axis(ax)

% and uncertainty
%--------------------------------------------------------------------------
subplot(2,1,2)

MDP.C = C;
c     = linspace(.5,.9,4);
for i = 1:length(c)
    
    % preference
    %----------------------------------------------------------------------
    a      = c(i);
    U{1,1} = [.5 .5; .5 .5; 0 0; 0 0];
    U{2,2} = [0 0; 0 0; a (1 - a); (1 - a) a];
    U{3,3} = [0 0; 0 0; (1 - a) a; a (1 - a)];
    U{4,4} = [1 0; 0 1; 0 0; 0 0];
    MDP.A  = spm_cat(U);
    
    % simulate trials and record outcomes
    %----------------------------------------------------------------------
    MDP  = spm_MDP_game(MDP,'FE');
    col  = [1 1 1]*(length(c) - i)/length(c);
    plot(MDP.da,'Color',col), hold on
    
end
title('Uncertainty','FontSize',16)
xlabel('Peristimulus time','FontSize',12)
ylabel('Response','FontSize',12)
axis square, axis(ax)
