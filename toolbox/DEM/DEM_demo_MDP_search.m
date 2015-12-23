function MDP = DEM_demo_MDP_search
% Demo of active inference for visula salience
%__________________________________________________________________________
%
% This routine uses a Markov decision process formulation of active
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
% (that take the agent to the four locations) and four exteroceptive 
% outcomes (that depend on the agents locations) plus three interoceptive
% outcomes indicating reward (or not).
%
% This version focuses on factorising the hidden states causing
% (factorised) outcomes. This factorisation is implicit in the tensor
% production used in the companion demo.  Here the factorisation is explicit
% enabling us to model multiple modalities (outcome factors) and distinct
% hidden causes of observation (hidden state factors like what and where).
% The behaviour is formally similar to the vanilla scheme but allows a much
% more intuitive (and possibly flexible) model specification.
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_MDP_search.m 6655 2015-12-23 20:21:27Z karl $
 
% set up and preliminaries
%==========================================================================
rng('default')

% prior beliefs about initial states (in terms of counts_: D and d
%--------------------------------------------------------------------------
d{1} = [8 8 8]';           % what:  {'flee','feed','wait'}
d{2} = [1 0 0 0 0 0 0 0]'; % where: {'start','1',...,'4','flee','feed','wait'}
d{3} = [8 8]';             % flip:  {'up','down'}
d{4} = [8 8]';             % flip:  {'left','right'}


% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(d);
for f = 1:Nf
    Ns(f) = numel(d{f});
end
No    = [6 8];
Ng    = numel(No);
for g = 1:Ng
    A{g} = zeros([No(g),Ns]);
end
for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
        for f3 = 1:Ns(3)
            for f4 = 1:Ns(4)
                
                % latent cues for this hidden state
                %----------------------------------------------------------
                if f1 == 1, a = {'bird','cat' ;'null','null'}; end
                if f1 == 2, a = {'bird','seed';'null','null'}; end
                if f1 == 3, a = {'bird','null';'null','seed'}; end
                
                if f3 == 2, a = flipud(a); end
                if f3 == 2, a = fliplr(a); end
                
                % what: A{1} {'null','bird,'seed','cat','right','wrong'}
                %----------------------------------------------------------
                if f2 == 1
                    A{1}(1,f1,f2,f3,f4) = 1;           
                elseif f2 > 1 && f2 < 5
                    A{1}(1,f1,f2,f3,f4) = strcmp(a{f2 - 1},'null');
                    A{1}(2,f1,f2,f3,f4) = strcmp(a{f2 - 1},'bird');
                    A{1}(3,f1,f2,f3,f4) = strcmp(a{f2 - 1},'seed');
                    A{1}(4,f1,f2,f3,f4) = strcmp(a{f2 - 1},'cat');
                else
                    A{1}(5,f1,f2,f3,f4) = (f2 - 5) == f1;
                    A{1}(6,f1,f2,f3,f4) = (f2 - 5) ~= f1;
                end
                
                % where: A{2} {'start','1',...,'4','flee','feed','wait'}
                %----------------------------------------------------------
                A{2}(f2,f1,f2,f3,f4) = 1;
                
            end
        end
    end
end
for g = 1:Ng
    A{g} = double(A{g});
end
 
% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% controlable fixation points (
%--------------------------------------------------------------------------
for k = 1:Ns(2)
    B{2}(:,1:5,k)   = 0;
    B{2}(k,1:5,k)   = 1;        % move to the k-th location
    B{2}(6:8,6:8,k) = eye(3);   % last three locations are absorbing
end


% allowable policies (of depth T): sequences of actions for each factor
%--------------------------------------------------------------------------
Nu       = 8;
U(1,:)   = kron(1:Nu,ones(1,Nu));
U(2,:)   = kron(ones(1,Nu),1:Nu);
T        = size(U,1);
Np       = size(U,2);
V        = ones(T,Np,Nf);
V(:,:,2) = U;


% priors: (utility) C
%--------------------------------------------------------------------------
T         = size(V,1);
c         = 4;
C{1}      = zeros(No(1),T);
C{2}      = zeros(No(2),T);
C{1}(5,:) =  c;                  % the agent expects to be right
C{1}(6,:) = -c;                  % and not wrong
 
 
% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.V = V;                    % allowable policies
mdp.A = A;                    % observation model
mdp.B = B;                    % transition probabilities
mdp.C = C;                    % preferred outcomes
mdp.d = d;                    % prior over initial states
mdp.s = ones(Nf,1);           % true initial state
 
mdp.Aname = {'what','where'};
mdp.Bname = {'what','where','flip','flip'};

MDP  = spm_MDP_VB_X(mdp);
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);

return
 
% true initial states – with context change at trial 12
%--------------------------------------------------------------------------
i              = [1,3];          % change context in a couple of trials
[MDP(1:32)]    = deal(mdp);      % create structure array
[MDP(i).s]     = deal([1 2]');   % deal context changes
 
 
% Solve - an example game: a run of reds then an oddball
%==========================================================================
MDP  = spm_MDP_VB_X(MDP);
 
% illustrate behavioural responses – first trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1));
 
% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_game(MDP);
 
% illustrate phase-precession and responses to chosen option - 1st trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_LFP(MDP(1),[2 3;3 3],1);
 
% illustrate phase-amplitude (theta-gamma) coupling
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MDP(1:8));
 
% illustrate familiarity (c.f., MMN) and context learning
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
spm_MDP_VB_LFP(MDP([2,16]),[1;2],2);
subplot(4,1,1), title('Repetition suppression and DA transfer','FontSize',16)
 
spm_figure('GetWin','Figure 6'); clf
v  = spm_MDP_VB_LFP(MDP([2,16]),[1;2],2);
t  = (1:16)*16 + 80;
subplot(2,1,1),plot(t,v{1}{2,1},'b-.',t,v{2}{2,1},'b:',t,v{2}{2,1} - v{1}{2,1})
xlabel('Time (ms)'),ylabel('LFP'),title('Difference waveform (MMN)','FontSize',16)
legend({'oddball','standard','MMN'}), grid on, axis square


