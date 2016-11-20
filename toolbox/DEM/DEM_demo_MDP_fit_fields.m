function DCM = DEM_demo_MDP_fit_fields
% Demo of active inference for visual salience
%__________________________________________________________________________
%
% This routine uses active inference for Markov decision processes to
% illustrate epistemic foraging in the context of visual searches. Here,
% the agent has to categorise scenes on the basis of the relative position
% of various cues. Crucially, the agent can only sample one cue or location
% at a time and therefore has to accumulate evidence for competing
% hypotheses. This rests upon resolving uncertainty about which scene or
% hypothesis is in play through the minimisation of expected free energy.
%
% When the agent become sufficiently confident about the underlying scene,
% it then makes a saccade to a choice location - to obtain feedback (right
% or wrong). The agent prefers to be right and does not expect to be
% wrong. We first illustrate a single trial in terms of behaviour and
% electrophysiological responses. We then consider sequences of trials and
% how one can recover prior preferences by inverting a model of observed
% responses (and cues).
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_MDP_fit_fields.m 6939 2016-11-20 13:33:56Z karl $



% create MDP structure of this paradigm 
%==========================================================================

% add prior preferences (and alpha) for this subject
%--------------------------------------------------------------------------
P.C      = log([2,1/4]);                 % log scaling for preferences
P.alpha  = log(16);                      % log scaling for alpha
mdp      = spm_MDP_gen(P);

% An example trial:
%--------------------------------------------------------------------------
SDP      = spm_MDP_VB_X(mdp);

% Illustrate belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(SDP);

% and phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(SDP,[],1);


% generate a sequence of trials
%==========================================================================
% use an efficient experimental design with all three categories presented
% with and without horizontal flipping (see spm_MDP_create below)
%--------------------------------------------------------------------------
s     = [];
for i = 1:3
    for j = 1:2
        k    = ones(4,1);
        k(1) = i;
        k(4) = j;
        s(:,end + 1) = k;
    end
end

for i = 1:size(s,2)
    MDP(i)   = mdp;      % create structure array
    MDP(i).s = s(:,i);   % context
end

% solve sequence to get stimuli {MDP.o} and choices {MDP.u}
%--------------------------------------------------------------------------
MDP  = spm_MDP_VB_X(MDP);



% Invert to recover parameters (preferences and precision)
%==========================================================================
% Computational phenotyping entails the estimation of subject specific
% priors (and other parameters) that best explain observed choices (under
% experimentally determined observations). This is done by using a Newton
% scheme and Variational Laplace, based on a log likelihood. The log
% likelihood of choices or actions is determined by the posterior
% distribution over choices. This log likelihood is returned by an
% auxiliary function: spm_MDP_L. The parameterisation of a subject's prior
% beliefs is specified with another auxiliary function that is configured
% for any particular experiment: spm_MDP_field. This function sets the
% requisite fields in the MDP scheme (usually MDP.C) as a function of
% user-specified parameters, whose priors are M.pE and M.pC. in the example
% spm_MDP_field subfunction below, we have used to parameters to scale
% preferences about outcomes and where to look, and a third parameter to
% set alpha - the precision with which an action is selected from posterior
% beliefs about action. By design, we expect the first parameter (about
% preferred outcomes to be about one and alpha to be about two. As there
% were no preferences about a speeded response used, during the generation
% of the simulated data, we expect the second parameter to be negative
% (because it is a log scaling parameter).


% response model priors (usually prior preferences in MDP)
%--------------------------------------------------------------------------
pE.C      = [0,0];                  % log scaling for preferences
pE.alpha  = 0;                      % log scaling for alpha
pC        = eye(spm_length(pE))/32; % shrinkage priors

% complete model specification
%--------------------------------------------------------------------------
DCM.M.G   = @spm_MDP_gen;           % parameterisation of MDP (see below)
DCM.M.L   = @spm_MDP_L;             % log-likelihood function
DCM.M.pE  = pE;                     % prior expectations
DCM.M.pC  = pC;                     % prior covariances
DCM.U     = {MDP.o};                % trial specification (stimuli)
DCM.Y     = {MDP.u};                % responses (action or choices)

% model inversion with Variational Laplace
%--------------------------------------------------------------------------
[Ep,Cp,F] = spm_nlsi_Newton(DCM.M,DCM.U,DCM.Y);

% Store posterior densities and log evidnce (free energy)
%--------------------------------------------------------------------------
DCM.Ep    = Ep;
DCM.Cp    = Cp;
DCM.F     = F;

% superimpose true values on posterior estimates
%--------------------------------------------------------------------------
p         = spm_vec(P);
subplot(2,2,1),hold on; plot(p(1),p(2),'.g','MarkerSize',32); hold off
subplot(2,2,4),hold on; plot(p,        '.g','MarkerSize',16); hold off

return


function MDP = spm_MDP_gen(pE)
% auxiliary function for creating MDP models
% FORMAT MDP = spm_MDP_gen(pE)
% pE  - subject specific  parameters
% MDP - MDP structure  for simulating choice behaviour
%
% This subfunction fills in the parameter fields of a MDP using free
% parameters specified in pE
%__________________________________________________________________________
%
% In this example, an agent has to categorise a scene that comprises
% potential cues at four peripheral locations, starting from a central
% fixation point. This involves a form of scene construction, in which the
% relationship between various cues determines the category of scene. In
% brief, the scene always contains a bird and seed, or bird and a cat. If
% the bird is next to the seed or the cat, then the scene is categorised as
% feed or flee respectively. Conversely, if the seed is in an opposite
% diagonal location, the category is wait. The particular positions of the
% cues are irrelevant, the important attribute are there relationships.
% This means hidden states have to include some spatial mappings that
% induce invariances to spatial transformations. These are reflections
% around the horizontal and vertical axes.
%
% There are two outcome modalities (what and where), encoding one of six
% cues and one of eight locations (there are three extra locations that
% provide feedback about the respective decision).  The hidden states
% have four factors; corresponding to context (the three categories), where
% (the eight locations) and two further factors modelling invariance
%--------------------------------------------------------------------------

% prior beliefs about initial states (in terms of counts_: D and d
%--------------------------------------------------------------------------
D{1} = [1 1 1]';           % what:     {'flee','feed','wait'}
D{2} = [1 0 0 0 0 0 0 0]'; % where:    {'start','1',...,'4','flee','feed','wait'}
D{3} = [1 1]';             % flip(ud): {'no','yes'}
D{4} = [1 1]';             % flip(rl): {'no','yes'}

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
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
                
                % location of cues for this hidden state
                %----------------------------------------------------------
                if f1 == 1, a = {'bird','cat' ;'null','null'}; end
                if f1 == 2, a = {'bird','seed';'null','null'}; end
                if f1 == 3, a = {'bird','null';'null','seed'}; end
                
                % flip cues according to hidden (invariants) states
                %----------------------------------------------------------
                if f3 == 2, a = flipud(a); end
                if f4 == 2, a = fliplr(a); end
                
                % what: A{1} {'null','bird,'seed','cat','right','wrong'}
                %==========================================================
                if f2 == 1
                    
                    % at fixation location
                    %----------------------------------------------------------
                    A{1}(1,f1,f2,f3,f4) = true;
                    
                elseif f2 > 1 && f2 < 6
                    
                    % saccade to cue location
                    %----------------------------------------------------------
                    A{1}(1,f1,f2,f3,f4) = strcmp(a{f2 - 1},'null');
                    A{1}(2,f1,f2,f3,f4) = strcmp(a{f2 - 1},'bird');
                    A{1}(3,f1,f2,f3,f4) = strcmp(a{f2 - 1},'seed');
                    A{1}(4,f1,f2,f3,f4) = strcmp(a{f2 - 1},'cat');
                    
                elseif f2 > 5
                    
                    % saccade choice location
                    %------------------------------------------------------
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

% controllable fixation points: move to the k-th location
%--------------------------------------------------------------------------
for k = 1:Ns(2)
    B{2}(:,:,k) = 0;
    B{2}(k,:,k) = 1;
end

% plus a bespoke state action (reading) policy
%--------------------------------------------------------------------------
B{2}(:,:,end + 1) = [ ...
    0     0     0     0     1     1     1     1
    1     0     0     0     0     0     0     0
    0     1     0     0     0     0     0     0
    0     0     1     0     0     0     0     0
    0     0     0     1     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0
    0     0     0     0     0     0     0     0];
 

% allowable policies (here, specified as the next action) U
%--------------------------------------------------------------------------
Np       = size(B{2},3);
U        = ones(1,Np,Nf);
U(:,:,2) = 1:Np;


% priors over policies (with a slight preference for reading)
%--------------------------------------------------------------------------
E        = ones(Np,1);
E(end)   = 1 + 1/8;


% priors: (utility) C
%--------------------------------------------------------------------------
T        = 5;
C{1}     = zeros(6,T);
C{2}     = zeros(8,T);

% patterns of preferences
%--------------------------------------------------------------------------
C{1}(5,:)   =  1;                 % the agent expects to be right
C{1}(6,:)   = -2;                 % and not wrong
C{2}(1:5,3) = -1;                 % make tardy sampling costly
C{2}(1:5,4) = -2;                 % make tardy sampling costly
C{2}(1:5,5) = -3;                 % make tardy sampling costly

% scale by parameters
%--------------------------------------------------------------------------
C{1}  = C{1}*exp(pE.C(1));
C{2}  = C{2}*exp(pE.C(2));


% MDP Structure
%==========================================================================
MDP.T = T;                      % number of moves
MDP.U = U;                      % allowable policies
MDP.A = A;                      % observation model
MDP.B = B;                      % transition probabilities
MDP.C = C;                      % preferred outcomes
MDP.D = D;                      % prior over initial states
MDP.E = E;                      % prior over initial states
MDP.s = [1 1 1 1]';             % initial state (flee)
MDP.o = [1 1]';                 % initial outcome

MDP.Aname = {'what','where'};
MDP.Bname = {'what','where','flip','flip'};

% alpha
%--------------------------------------------------------------------------
MDP.alpha = exp(pE.alpha);
