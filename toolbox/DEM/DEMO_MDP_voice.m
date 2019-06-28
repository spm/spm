function MDP = DEMO_MDP_voice
% Active inference in conversation
%__________________________________________________________________________
%
% This routine provide simulations of reading to demonstrate deep temporal
% generative models. It builds upon the scene construction simulations to
% equip the generative model with a second hierarchical level. In effect,
% this creates an agent that can accumulate evidence at the second level
% based upon epistemic foraging at the first. In brief, the agent has to
% categorise a sentence or narrative into one of two categories (happy or
% sad), where it entertains six possible sentences. Each sentence comprises
% four words, which are themselves constituted by two pictures or graphemes
% These are the same visual outcomes used in previous illustrations of
% scene construction and saccadic searches.
%
% Here, the agent has policies at two levels. The second level policy (with
% just one step into the future) allows it to either look at the next word
% or stay on the current page and make a decision. Concurrently, a first
% level policy entails one of four saccadic eye movements to each quadrant
% of the current page, where it will sample a particular grapheme.
%
% This provides a rough simulation of reading – that can be made more
% realistic by terminating first level active inference, when there can be
% no further increase in expected free energy (i.e., all uncertainty about
% the current word has been resolved). The subsequent inferred hidden
% states then become the outcome for the level above.
%
% To illustrate the schemes biological plausibility, one can change the
% agent’s prior beliefs and repeat the reading sequence under violations of
% either local (whether the graphemes are flipped vertically) or globally
% (whether the sentence is surprising) expectations. This produces a
% mismatch negativity (MMN) under local violations) and a MMN with a
% P300 with global violations.
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEMO_MDP_voice.m 7631 2019-06-28 12:25:34Z karl $
 
 
% rng('default')
 
% first level (lexical)
%==========================================================================
% question 1: {'noun (7)'}
% question 2: {'noun (7)','adverb (9)'}
% question 3: {'adjective (8)','noun (7)','adverb (9)'}
%--------------------------------------------------------------------------
% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------

% setup lexicon and syntactic (sentence) structures
%==========================================================================
% note each sequence ends in the same (space) outcome: s010 = ' '
syn{1}      = {'is'};
syn{2}      = {'there'};
syn{3}      = {'a'};
syn{4}      = {'yes','OK'};
syn{5}      = {'no','sorry'};
syn{6}      = {'Iam'};
syn{7}      = {'not'};
syn{8}      = {'sure'};
syn{9}      = {'ready','OK','go'};
syn{10}     = {' '};

noun        = {'square','triangle'};
adjective   = {'green','red'};
adverb      = {'above','below'};

% allowable outcomes (words)
%--------------------------------------------------------------------------
outcome     = unique([syn{:} noun{:} adjective{:} adverb{:}]);

% allowable sentences
%--------------------------------------------------------------------------
sentence{1} = {'s11','s12','s13','&1f1','s010'};
sentence{2} = {'s21','s22','s23','&2f1','&2f3','s010'};
sentence{3} = {'s31' 's33','&3f2','&3f1','&3f3','s010'};
sentence{4} = {'s44','s010'};
sentence{5} = {'s55','s010'};
sentence{6} = {'s66','s67','s68','s010'};
sentence{7} = {'s79','s010'};

% assemble hidden (syntax) states
%--------------------------------------------------------------------------
label.factor{1} = 'noun';   label.name{1}  = noun;
label.factor{2} = 'adj.';   label.name{2}  = adjective;
label.factor{3} = 'adverb'; label.name{3}  = adverb;
label.factor{4} = 'syntax'; label.name{4}  = unique([sentence{:}]);

% prior beliefs about initial states D 
%--------------------------------------------------------------------------
for i = 1:numel(label.factor)
    n    = numel(label.name{i});
    D{i} = ones(n,1)/n;
end

% restrict initial states to the beginning of a sentence
%--------------------------------------------------------------------------
state = label.name{4};
D{4}  = spm_zeros(D{4});
for i = 1:numel(sentence)
    j = ismember(state,sentence{i}(1));
    D{4}(j) = 1;
end

% hidden factors
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end

% single outcome modality with multiple phrases
%--------------------------------------------------------------------------
label.modality{1} = 'word';   label.outcome{1} = outcome;

for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
        for f3 = 1:Ns(3)
            for f4 = 1:Ns(4)
                
                % indices of the state
                %----------------------------------------------------------
                j = {f1,f2,f3,f4};
                
                % words under this state
                %==========================================================
                name = label.name{4}{f4};
                if name(1) == '&'
                    
                    % get noun, adverb or adjective
                    %------------------------------------------------------
                    p   = eval(name(end));
                    q   = eval(name(3:end));
                    out = label.name{p}(q);
                    
                elseif name(1) == 's'
                    
                    % otherwise get synonyms
                    %------------------------------------------------------
                    q   = eval(name(3:end));
                    out = syn{q};
                end
                
                % place words in likelihood matrix
                %----------------------------------------------------------
                i = ismember(outcome,out);
                A{1}(i,j{:}) = 1/sum(i);
                
                % prosidy under this state
                %==========================================================
                for p = 2:4
                    A{p}(:,j{:}) = full(sparse(3:6,1,1,8,1));
                end

            end
        end
    end
end

% transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end
 
% specify syntax; i.e., the sequence of syntactic (sentence) states
%--------------------------------------------------------------------------
B{4}  = spm_zeros(B{4});
for s = 1:numel(sentence)
    for t = 2:numel(sentence{s})
        i = find(ismember(state,sentence{s}(t - 1)));
        j = find(ismember(state,sentence{s}(t)));
        B{4}(j,i) = 1;
    end
    B{4}(j,j) = 1;
end

 
% MDP Structure
%--------------------------------------------------------------------------
mdp.T = 6;                      % number of updates
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.D = D;                      % prior over initial states
mdp.o = [];

mdp.VOX   = [0,0,0];
mdp.label = label;
MDP       = spm_MDP_check(mdp);

clear A B D mdp


%% second level (narrative)
%==========================================================================
% question 1: {'noun (7)'}
% question 2: {'noun (7)','adverb (9)'}
% question 3: {'adjective (8)','noun (7)','adverb (9)'}
%--------------------------------------------------------------------------
label.factor{1}  = 'narrative';    label.name{1}  = {'ready','question','answer'};
label.factor{2}  = 'question';     label.name{2}  = {'is?','where?','what?'};
label.factor{3}  = 'upper colour'; label.name{3}  = {'green','red'};
label.factor{4}  = 'lower colour'; label.name{4}  = {'green','red'};
label.factor{5}  = 'upper shape';  label.name{5}  = {'square','triangle'};
label.factor{6}  = 'lower shape';  label.name{6}  = {'square','triangle'};
label.factor{7}  = 'noun';         label.name{7}  = {'square','triangle'};
label.factor{8}  = 'adjective';    label.name{8}  = {'green','red'};
label.factor{9}  = 'adverb';       label.name{9}  = {'above','below'};

% prior beliefs about initial states D 
%--------------------------------------------------------------------------
for i = 1:numel(label.factor)
    n    = numel(label.name{i});
    D{i} = ones(n,1)/n;
end

% known initial states
%--------------------------------------------------------------------------
D{1}(1)  = 128;

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
label.modality{1} = 'noun';   label.outcome{1}  = MDP.label.name{1};
label.modality{2} = 'adj.';   label.outcome{2}  = MDP.label.name{2};
label.modality{3} = 'adverb'; label.outcome{3}  = MDP.label.name{3};
label.modality{4} = 'syntax'; label.outcome{4}  = MDP.label.name{4}(find(MDP.D{4}));

Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
for f1 = 1:Ns(1) % narrative {'ready','question','answer'}
    for f2 = 1:Ns(2) % question {'is?','where?','what?'}
        for f3 = 1:Ns(3) % upper colour {'green','red'}
            for f4 = 1:Ns(4) % lower colour {'green','red'}
                for f5 = 1:Ns(5) % upper shape {'square','triangle'}
                    for f6 = 1:Ns(6) % lower shape {'square','triangle'}
                        for f7 = 1:Ns(7) % noun {'square','triangle'}
                            for f8 = 1:Ns(8) % adjective {'green','red'}
                                for f9 = 1:Ns(9) % adverb {'above','below'}
                                    
                                    % indices
                                    %--------------------------------------
                                    j = {f1,f2,f3,f4,f5,f6,f7,f8,f9};
                                    
                                    % answer: depending on question and beliefs
                                    %--------------------------------------
                                    Y  = 0;
                                    if f1 == 3
                                        if f2 == 1
                                            Y = (f7 == f5) | (f7 == f6);
                                        elseif f2 == 2
                                            if f9 == 1
                                                Y = (f7 == f5);
                                            elseif f9 == 2
                                                Y = (f7 == f6);
                                            end
                                        elseif f2 == 3
                                            if f9 == 1
                                                Y = (f7 == f5) & (f8 == f3);
                                            elseif f9 == 2
                                                Y = (f7 == f6) & (f8 == f4);
                                            end
                                        end
                                    end
                                    
                                    % A{1} noun:
                                    %======================================
                                    A{1}(f7,j{:}) = 1;
                                    
                                    % A{2} adjective:
                                    %======================================
                                    A{2}(f8,j{:}) = 1;
                                    
                                    % A{3} adverb:
                                    %======================================
                                    A{3}(f9,j{:}) = 1;
                                    
                                    % A{4} syntax: {Q1,Q2,Q3,Y,N,?,Go}
                                    %======================================
                                    if f1 == 1
                                        A{4}(7,j{:})     = 1;
                                    elseif f1 == 2
                                        A{4}(f2,j{:})    = 1;
                                    elseif f1 == 3
                                        if Y
                                            A{4}(4,j{:}) = .98;
                                        else
                                            A{4}(5,j{:}) = .98;
                                        end
                                        A{4}(6,j{:})     = .02;
                                    end
                                    
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

% likelihood mappings
%--------------------------------------------------------------------------
Ng    = numel(A);
for g = 1:Ng
    No(g) = size(A{g},1);
end
 
% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
    label.action{f} = {'stay'};
end
 
% transitions B(1): {'ready','question','answer'}
%--------------------------------------------------------------------------
B{1}(:,:,1) = spm_speye(Ns(1),Ns(1),-1); B{1}(1,Ns(1),1) = 1;

% control states B(2): question {'1','2' or '3'} & control states B(7):B{9} 
%--------------------------------------------------------------------------
% D{2} ;      % question;     {'1','2','3'};
% D{7} ;      % noun:         {'square','triangle'}
% D{8} ;      % adjective:    {'green','red'}
% D{9} ;      % adverb:       {'above','below'}
%--------------------------------------------------------------------------
for f = [2 7 8 9]
    for k = 1:Ns(f)
        B{f}(:,:,k + 1) = 0;
        B{f}(k,:,k + 1) = 1;
    end
    label.action{f} = ['stay', label.name{f}];
end

% allowable policies (time x polcy x factor): ready, question and answer
%--------------------------------------------------------------------------
% question 1: {'noun (7)'}
% question 2: {'noun (7)','adverb (9)'}
% question 3: {'adjective (8)','noun (7)','adverb (9)'}
%--------------------------------------------------------------------------
V         = ones(2,14,Nf);
V(1,:,2)  = 1 + [1 1 2 2 2 2 3 3 3 3 3 3 3 3];
V(1,:,7)  = 1 + [1 2 1 1 2 2 1 1 1 1 2 2 2 2];
V(1,:,8)  = 1 + [1 1 1 1 1 1 1 1 2 2 1 1 2 2];
V(1,:,9)  = 1 + [1 1 1 2 1 2 1 2 1 2 1 2 1 2];
V(2,:,:)  = 1;

 
% priors: (utility) C: A{4} syntax: {'1','2','3','Y','N','?'}
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end
C{4}(4,:) =  1/4;              % and affirmative answers
C{4}(5,:) = -1/4;              % and negative answers
 
% actual state of the world
%--------------------------------------------------------------------------
s    = ones(Nf,1);
s(4) = 2;


% MDP Structure
%--------------------------------------------------------------------------
mdp.FCN    = @spm_questions_plot;
mdp.MDP    = MDP;
mdp.label  = label;             % names of factors and outcomes
mdp.tau    = 4;                 % time constant of belief updating
mdp.erp    = 4;                 % initialization

mdp.T = 3;                      % ready, question, answer
mdp.V = V;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states (context)
mdp.s = s;                      % initial state
mdp.o = [];                     % outcomes

mdp.link = spm_MDP_link(mdp);   % map outputs to initial (lower) states

 
%% illustrate questioning
%==========================================================================
clear MDP
[mdp.MDP(1,1:3)] = deal(mdp.MDP);

% agent asks (by setting VOX to [0 1 0]
%--------------------------------------------------------------------------
VOX   = [2, 1, 0];
for t = 1:mdp.T
    mdp.MDP(t).VOX = VOX(t);
end
OPTIONS.D    = 1;
[MDP(1,1:6)] = deal(mdp);


% agent asks (by setting outomes to [0 1 0]
%--------------------------------------------------------------------------

% belief updating
%==========================================================================
spm_figure('GetWin','20 Questions'); clf; subplot(4,3,2); axis off
MDP   = spm_MDP_VB_X(MDP,OPTIONS);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1),[1 2 4],[1 3 4]);

% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],4);

spm_figure('GetWin','Figure 2A'); clf
spm_MDP_VB_LFP(MDP,[],4,1);

spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_ERP(MDP(4:6),[4,3]);

spm_figure('GetWin','20 Questions'); clf
for i = 1:size(MDP,2)
    subplot(4,3,i)
    spm_questions_plot(MDP(1,i))
end

return



function spm_questions_plot(MDP,X)
%% illustrate beliefs
%--------------------------------------------------------------------------

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
% label.modality{1} = 'noun';   label.outcome{1}  = {'square','triangle'};
% label.modality{2} = 'adj.';   label.outcome{2}  = {'green','red'};
% label.modality{3} = 'adverb'; label.outcome{3}  = {'above','below'};
% label.modality{4} = 'syntax'; label.outcome{4}  = {'1','2','3','Y','N','Sil.'};

% plot question, answer and posterior belief
%==========================================================================
spm_figure('GetWin','20 Questions'); hold off
if nargin < 2; X = MDP.X; end

% Assemble question-and-answer
%----------------------------------------------------------------------
answer    = MDP.o(4,3);
if ~answer
    return
end

qstr      = MDP.MDP(1).label.outcome{1}(MDP.mdp(2).o(1,:));
qstr(2,:) = {' '}; qstr = cell2mat(qstr(:)');
astr      = MDP.MDP(1).label.outcome{1}(MDP.mdp(3).o(1,:));
astr(2,:) = {' '}; astr = cell2mat(astr(:)');


% is the answer right (for a single player)?
%----------------------------------------------------------------------
ind    = num2cell(MDP.s(:,3));
if answer == find(MDP.A{4}(:,ind{:}),1)
    cor = spm_softmax(2*[0;1;0]);
else
    cor = spm_softmax(2*[1;0;0]);
end

%  plot posterior beliefs
%----------------------------------------------------------------------
% label.factor{3}  = 'upper colour'; label.name{3}  = {'green','red'};
% label.factor{4}  = 'lower colour'; label.name{4}  = {'green','red'};
% label.factor{5}  = 'upper shape';  label.name{5}  = {'square','triangle'};
% label.factor{6}  = 'lower shape';  label.name{6}  = {'square','triangle'};

% upper and lower object
%----------------------------------------------------------------------
T       = 2;
col{1}  = [X{3}(2,T) X{3}(1,T) 0];
col{1}  = col{1}*X{5}(1,T) + (1 - X{5}(1,T));
col{2}  = [X{3}(2,T) X{3}(1,T) 0];
col{2}  = col{2}*X{5}(2,T) + (1 - X{5}(2,T));
col{3}  = [X{4}(2,T) X{4}(1,T) 0];
col{3}  = col{3}*X{6}(1,T) + (1 - X{6}(1,T));
col{4}  = [X{4}(2,T) X{4}(1,T) 0];
col{4}  = col{4}*X{6}(2,T) + (1 - X{6}(2,T));

plot(1,0,'^','MarkerSize',24,'LineWidth',4,'Color',col{4}), hold on
plot(1,0,'s','MarkerSize',24,'LineWidth',4,'Color',col{3})
plot(1,1,'^','MarkerSize',24,'LineWidth',4,'Color',col{2})
plot(1,1,'s','MarkerSize',24,'LineWidth',4,'Color',col{1})

text(1, 2,qstr,'HorizontalAlignment','Center')
text(1,-1,astr,'HorizontalAlignment','Center','FontWeight','bold','Color',cor)
axis([0 2 -1.5 2.5]), axis off, axis square

% upper and lower object
%--------------------------------------------------------------------------
rgb     = {[0 1 0],[1 0 0]};
shape   = {'s','^'};
plot(1 + 1/2,1,shape{MDP.s(5)},'MarkerSize',8,'LineWidth',1,'Color',rgb{MDP.s(3)})
plot(1 + 1/2,0,shape{MDP.s(6)},'MarkerSize',8,'LineWidth',1,'Color',rgb{MDP.s(4)})
drawnow





