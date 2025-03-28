function MDP = DEM_laws
% Demo of context sensitive preferences
%__________________________________________________________________________
% 
% This demonstration script illustrates context-sensitive preferences and
% the resolution of conflicting preferences under active inference. It
% illustrates how one set of preferences (e.g., preferring not to violate
% social or institutional norms) can be overwritten by potentially
% conflicting preferences (e.g., preferred outcomes in an emergency), in
% terms of choice behaviour. The particular setup considered here involves
% moving from a starting location (e.g., the sidewalk) to a street crossing
% and crossing the road, when a (green) walk sign allows. Crucially, a
% siren can sound at any time, indicating an imperative to reach the other
% side of the road. In nominal circumstances the agent will wait patiently
% until the walk sign (i.e., deontic cue) switches from don't walk (red) to
% walk (green). However, if the siren sounds, the agent could — with the
% right preferences – cross the road, risking aversive outcomes; such as
% being arrested or run over by a car.
%
% In the example below, the agent has mild preferences for being at her
% destination, which are not sufficient to overcome her aversion to
% violations, unless she hears a siren. The don't walk sign (that can only
% be seen from the crossing) is on for six timesteps. After four timesteps
% the siren sounds and the agent infers there is an emergency. This compels
% her to cross the road, thereby eliciting a violation but enabling her to
% reach the safe destination. Interestingly, when the agent is sufficiently
% confident there is an emergency — given an ambiguous (auditory)
% likelihood mapping — she resolves her uncertainty about what to do,
% eliciting a simulated phasic release of dopamine. This be regarded as
% reporting the affective aspect of belief updating. Technically, dopamine
% scores the increases or decreases in the average expected free energy
% over policies; intuitively, this corresponds to the confidence in what
% the agent — believes she — is doing.
%__________________________________________________________________________
 
% set up and preliminaries
%==========================================================================
rng('default')

% latent states and outcomes
%==========================================================================
label.factor   = {'Location', 'Law' ,'Emergency'};
label.name     = {{'1','2','3','4'}, {'stay','go'}, {' ','urgent'}};
label.modality = {'Location', 'Walk sign', 'Alarm'};
label.outcome  = {{'1','2','3','4'}, {'red','green'}, {' ','alarm'}, {' ','alert'}};
label.action   = {{'stay','go'}, {' '}, {' '}};

% size of latent states and outcomes
%==========================================================================

% states
%--------------------------------------------------------------------------
Ns(1) = 4;                                % location
Ns(2) = 2;                                % stay/go context
Ns(3) = 2;                                % nominal/urgent context

% outcomes
%--------------------------------------------------------------------------
No(1) = 4;                                % locations (e.g., GPS)
No(2) = 2;                                % stay/go cue (e.g., crossing sign)
No(3) = 2;                                % emergency cue (e.g., siren)
No(4) = 2;                                % violation cue (e.g., car horn)

% action
%--------------------------------------------------------------------------
Nu = [2,1,1];                             % location is controllable

% initialise likelihoods A and priors B
%==========================================================================
Nf = numel(Ns);
Ng = numel(No);
for g = 1:Ng
    A{g}    = zeros([No(g),Ns]);
    C{g}    = zeros([No(g),Ns]);
    id.A{g} = 1:Nf;
    id.C{g} = 1:Nf;
end
for f = 1:Nf
    B{f} = zeros(Ns(f),Ns(f),Nu(f));
    D{f} = zeros(Ns(f),1);
end

% specify generative model under combinations of hidden states
%--------------------------------------------------------------------------
a     = 1/8;                              % a small uncertainty parameter
S     = spm_combinations(Ns);
for i = 1:size(S,1)

    % hidden states
    %----------------------------------------------------------------------
    loc  = S(i,1);                        % location
    sta  = S(i,2);                        % stay/go
    nom  = S(i,3);                        % nominal/urgent

    % likelihoods and preferences (A and C)
    %======================================================================

    % location (precise mapping)
    %----------------------------------------------------------------------
    A{1}(loc,loc,sta,nom) = 1;

    % slight preference for destination that increases under urgency
    %----------------------------------------------------------------------
    if nom == 1
        C{1}(:,loc,sta,nom) = spm_softmax([0;0;0;1/4]);
    else
        C{1}(:,loc,sta,nom) = spm_softmax([0;0;0;3]);
    end

    % stay/go cue (only seen from locations 2 and 3)
    %----------------------------------------------------------------------
    if loc == 2 || loc == 3
        A{2}(sta,loc,sta,nom) = 1;
    else
        A{2}(:,loc,sta,nom) = [1; 1]/2;
    end

    % emergency cue (heard everywhere, with slight ambiguity)
    %----------------------------------------------------------------------
    if nom == 1
        A{3}(:,loc,sta,nom) = [(1 - a);a];
    else
        A{3}(:,loc,sta,nom) = [a;(1 - a)];
    end

    % violation (when crossing under a stay context)
    %----------------------------------------------------------------------
    if loc == 3 && sta == 1
        A{4}(2,loc,sta,nom) = 1;
    else
        A{4}(1,loc,sta,nom) = 1;
    end

    % aversion for violation cue
    %----------------------------------------------------------------------
    C{4}(:,loc,sta,nom) = spm_softmax([3;0]);


    % Priors and initial conditions (B and D)
    %======================================================================
    
    % location (move to next location if u = 2, otherwise stay)
    %----------------------------------------------------------------------
    nex = loc + 1; if nex > Ns(1), nex = 3; end
    B{1}(loc,loc,1) = 1;
    B{1}(nex,loc,2) = 1;
    D{1}(1) = 1;                      % start at first location

    % stay/go (context could change at any time)
    %----------------------------------------------------------------------
    if sta == 1
        B{2}(:,sta) = [(1 - a);a];
    else
        B{2}(:,sta) = [a;(1 - a)];
    end
    D{2}   = [1; 1]/2;

    % Nominal/urgent context (could start at any time)
    %----------------------------------------------------------------------
    if nom == 1
        B{3}(:,nom) = [(1 - a);a];
    else
        B{3}(:,nom) = [0; 1];
    end
    D{3}   = [(1 - a); a];

end

% Create MDP structure
%==========================================================================
U     = Nu > 1;                       % controllable factors
MDP.T = 10;                           % number of moves
MDP.U = U;                            % controllable factors
MDP.A = A;                            % likelihood probabilities
MDP.B = B;                            % transition probabilities
MDP.C = C;                            % prior preferences
MDP.D = D;                            % prior over initial states
MDP.N = 3;                            % planning depth (4)

MDP.id    = id;                       % edges
MDP.label = label;                    % names

OPTIONS.N = 1;                        % switch on neuronal simulator

% Illustrate with an example by prescribing states of affairs
%==========================================================================
MDP.s      = zeros(Nf,MDP.T);
MDP.s(2,:) = [1 1 1 1 1 1 2 2 2 2];        % red light
MDP.s(3,:) = [1 1 1 2 2 2 2 2 2 2];        % siren

PDP = spm_MDP_VB_XXX(MDP,OPTIONS);
 
% illustrate active inference
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(PDP);

return
 
% illustrate behavioural responses and neuronal correlates
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_game(PDP);
 
% illustrate phase-precession and responses to chosen option - 1st trial
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_LFP(PDP);

return
