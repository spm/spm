function MDP = DEMO_MDP_questions
% Demo of active inference for visual salience
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
% This provides a rough simulation of reading - that can be made more
% realistic by terminating first level active inference, when there can be
% no further increase in expected free energy (i.e., all uncertainty about
% the current word has been resolved). The subsequent inferred hidden
% states then become the outcome for the level above.
%
% To illustrate the schemes biological plausibility, one can change the
% agent's prior beliefs and repeat the reading sequence under violations of
% either local (whether the graphemes are flipped vertically) or globally
% (whether the sentence is surprising) expectations. This produces a
% mismatch negativity (MMN) under local violations) and a MMN with a
% P300 with global violations.
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEMO_MDP_questions.m 7766 2020-01-05 21:37:39Z karl $
 
% set up and preliminaries: first level
%==========================================================================
%
%--------------------------------------------------------------------------
 
rng('default')
 
% first level (lexical)
%==========================================================================
% question 1: {'noun (7)'}
% question 2: {'noun (7)','adverb (9)'}
% question 3: {'adjective (8)','noun (7)','adverb (9)'}
%--------------------------------------------------------------------------
% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
sentence{1} = {'Is there a ','&1_1','?'};
sentence{2} = {'Is any ','&2_1','&2_3','?'};
sentence{3} = {'Is a ','&3_2','&3_1','&3_3','?'};
sentence{4} = {'Yes','!'};
sentence{5} = {'No','!'};
sentence{6} = {'Not sure','!'};
sentence{7} = {'Ready','?'};

syn{1}      = {'Ready','OK','please begin'};
syn{2}      = {'Not sure','I dont know','Cannot say'};
syn{3}      = {'No','Sorry, no','There is not'};
syn{4}      = {'Yes','Yes there is','Well done'};

% assemble hidden states and transition probabilities
%--------------------------------------------------------------------------
name  = {};
for i = 1:numel(sentence)
    for j = 1:numel(sentence{i})
        phrase = sentence{i}(j);
        if ~any(ismember(name,phrase))
        name(end + 1) = phrase;
        end
    end
end

label.factor{1} = 'noun';   label.name{1}  = {'square ','triangle '};
label.factor{2} = 'adj.';   label.name{2}  = {'green ','red '};
label.factor{3} = 'adverb'; label.name{3}  = {'above','below'};
label.factor{4} = 'syntax'; label.name{4}  = name;

% prior beliefs about initial states D 
%--------------------------------------------------------------------------
for i = 1:numel(label.factor)
    n    = numel(label.name{i});
    D{i} = ones(n,1)/n;
end

% restricted initial states to the beginning of a sentence
%--------------------------------------------------------------------------
state = label.name{4};
D{4}  = spm_zeros(D{4});
for i = 1:numel(sentence)
    j = ismember(state,sentence{i}(1));
    D{4}(j) = 1;
end

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end

outcome = {};
for i = 1:numel(label.name)
    outcome = [outcome,label.name{i}];
end
for i = 1:numel(syn)
    outcome = [outcome,syn{i}];
end
outcome = unique(outcome);
j     = [];
for i = 1:numel(outcome)
    if outcome{i}(1) ~= '&'
        j = [j,i];
    end
end
outcome = outcome(j);

% single outcome modality with multiple phrases
%--------------------------------------------------------------------------
label.modality{1} = 'phrase';   label.outcome{1} = outcome;

for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
        for f3 = 1:Ns(3)
            for f4 = 1:Ns(4)
                
                % indices of the state
                %----------------------------------------------------------
                j = {f1,f2,f3,f4};
                
                % assemble phrases under this state
                %----------------------------------------------------------
                name = label.name{4}{f4};
                if name(1) == '&'
                    p   = eval(name(end));
                    q   = eval(['f' name(end)]);
                    out = label.name{p}(q);
                else
                    out = name;
                    for i = 1:numel(syn)
                        if any(ismember(syn{i},name))
                            out = syn{i};
                        end
                    end
                end
                
                % placing likelihood matrix
                %----------------------------------------------------------
                i = ismember(outcome,out);
                A{1}(i,j{:}) = 1/sum(i);

            end
        end
    end
end

% transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end
 
% specifies syntax
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
mdp.T = 5;                      % number of updates
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.D = D;                      % prior over initial states
 
mdp.label = label;
mdp.chi   = 1/32;
 
clear A B D
 
MDP   = spm_MDP_check(mdp);

% % check belief updates (and behaviour)
% %------------------------------------------------------------------------
% MDP.s = [1 2 1 8]';
% MDP   = spm_MDP_VB_X(MDP);
%
% spm_figure('GetWin','Figure 1'); clf
% spm_MDP_VB_trial(MDP,[3 4],1);
% 
% % illustrate phase-precession and responses
% %------------------------------------------------------------------------
% spm_figure('GetWin','Figure 2'); clf
% spm_MDP_VB_LFP(MDP,[],4);
% 
% cell2mat(MDP.label.outcome{1}(MDP.o))


% set up and preliminaries: first level
%==========================================================================
% 
%--------------------------------------------------------------------------
 
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
for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
        for f3 = 1:Ns(3)
            for f4 = 1:Ns(4)
                for f5 = 1:Ns(5)
                    for f6 = 1:Ns(6)
                        for f7 = 1:Ns(7)
                            for f8 = 1:Ns(8)
                                for f9 = 1:Ns(9)
                                    
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
                                    
                                    % A{4} syntax: {'1','2','3','Y','N','?'}
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
                                         A{4}(6,j{:})    = .02;
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
    if any(spm_vec(sum(A{g})) ~= 1), disp('check A'), return, end
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

% allowable policies (time x polcy x factor)
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


% priors: (utility) C: A{4} syntax: {'1','2','3','Y','N','S'}
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end

% the agent expects affirmative answers
%--------------------------------------------------------------------------
C{4}(4,:) =  1/4;
C{4}(5,:) = -1/4;

% actual state of the world
%--------------------------------------------------------------------------
s    = ones(Nf,1);
s(4) = 2;

% MDP Structure
%--------------------------------------------------------------------------
mdp.MDP    = MDP;
mdp.label  = label;             % names of factors and outcomes
mdp.tau    = 4;                 % time constant of belief updating
mdp.erp    = 4;                 % initialization
mdp.chi    = 0;                 % initialization

mdp.V = V;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states (context)
mdp.s = s;                      % initial state
mdp.o = [];                     % outcomes

mdp.link = spm_MDP_link(mdp);
MDP      = spm_MDP_check(mdp);
%   spm_MDP_factor_graph(mdp);
 
 
%% illustrate questioning
%==========================================================================
clear MDP
OPTIONS.D    = 1;
[MDP(1,1:6)] = deal(mdp);

MDP   = spm_MDP_VB_X(MDP,OPTIONS);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP(1),[1 2 4],[1 3 4]);

% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],3);

spm_figure('GetWin','Figure 2A'); clf
spm_MDP_VB_LFP(MDP,[],3,1);

spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_ERP(MDP(4:6),[3,2]);

spm_figure('GetWin','Figure 4'); clf
for i = 1:size(MDP,2)
    subplot(4,3,i)
    spm_questions_plot(MDP(1,i))
end

% illustrate violations
%==========================================================================
j      = 5;                               % which answer
NDP    = MDP(j);                          % get states and outcomes
if NDP.o(4,3) == 4                        % switch the answer
    NDP.o(4,3) = 5;
else
    NDP.o(4,3) = 4;
end
NDP = rmfield(NDP,'link');                % remove link (outomes are given)
NDP = spm_MDP_VB_X(NDP);

% find greatest effect on belief updating (about the scene)
%--------------------------------------------------------------------------
for f = 3:6
    v(f) =  norm(spm_vec(MDP(j).X{f}) - spm_vec(NDP.X{f}),'inf');
end
[v,f] = max(v);

% responses to appropriate and inappropriate answers
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 6 expected' ); clf, spm_MDP_VB_LFP(MDP(5),[],f);
spm_figure('GetWin','Figure 7 violation'); clf, spm_MDP_VB_LFP(NDP   ,[],f);


% illustrate 'communication'
%==========================================================================

% increase efficiency (i.e., suppress neurophysiological correlates)
%--------------------------------------------------------------------------
mdp.tau    = 3;                         % time constant of belief updating
mdp.erp    = 1;                         % initialization
mdp.chi    = 1/64;                      % initialization

% create two agents
%--------------------------------------------------------------------------
clear MDP
[MDP(1:2,1:6)] = deal(mdp);

% give a second subject veridical beliefs about the scene
%--------------------------------------------------------------------------
for i = [3 4 5 6]
    MDP(2,1).D{i} = sparse(s(i),1,1,Ns(i),1);
end

% confident player answers then asks
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 8'); clf
for i = 1:length(MDP)
    
    if i < 5
        % first model asks and the second answers
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[2 1 2];
        MDP(2,i).n = ones(Ng,1)*[2 1 2];
    else
        % switch roles
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[1 2 1];
        MDP(2,i).n = ones(Ng,1)*[1 2 1];
    end
end

TDP   = spm_MDP_VB_X(MDP,OPTIONS);
for i = 1:size(TDP,2)
    subplot(4,3,i), spm_questions_plot(TDP(:,i))
end

% confident player asks then answers
%--------------------------------------------------------------------------
for i = 1:length(MDP)
    
    if i < 5
        % first model asks and the second answers
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[1 2 1];
        MDP(2,i).n = ones(Ng,1)*[1 2 1];
    else
        % switch roles
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[2 1 2];
        MDP(2,i).n = ones(Ng,1)*[2 1 2];
    end
end

TDP   = spm_MDP_VB_X(MDP,OPTIONS);
for i = 1:size(TDP,2)
    subplot(4,3,i + 6)
    spm_questions_plot(TDP(:,i))
end


% confident player tells then asks
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 9'); clf
for i = 1:length(MDP)
    
    if i < 7
        % first model asks and the second answers
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[1 2 2];
        MDP(2,i).n = ones(Ng,1)*[1 2 2];
    else
        % switch roles
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[1 2 1];
        MDP(2,i).n = ones(Ng,1)*[1 2 1];
    end
end

TDP   = spm_MDP_VB_X(MDP,OPTIONS);
for i = 1:size(TDP,2)
    subplot(4,3,i)
    spm_questions_plot(TDP(:,i))
end

% building fantasies
%==========================================================================
spm_figure('GetWin','Figure 10'); clf

% make both players uncertain about the scene at hand
%--------------------------------------------------------------------------
MDP(2,1).D = MDP(1,1).D;

for i = 1:length(MDP)
    
    if i < 3
        % first model asks and the second answers
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[1 2 1];
        MDP(2,i).n = ones(Ng,1)*[1 2 1];
    else
        % switch roles
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[2 1 2];
        MDP(2,i).n = ones(Ng,1)*[2 1 2];
    end
end

TDP   = spm_MDP_VB_X(MDP,OPTIONS);
for i = 1:size(TDP,2)
    subplot(4,3,i)
    spm_questions_plot(TDP(:,i))
end

% now repeat when precluding uncertain responses
%--------------------------------------------------------------------------
for i = 1:length(MDP)
    
    MDP(1,i).A{4} = MDP(1,i).A{4} > 1/2;
    MDP(2,i).A{4} = MDP(2,i).A{4} > 1/2;

    if i < 3
        % first model asks and the second answers
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[1 2 1];
        MDP(2,i).n = ones(Ng,1)*[1 2 1];
    else
        % switch roles
        %------------------------------------------------------------------
        MDP(1,i).n = ones(Ng,1)*[2 1 2];
        MDP(2,i).n = ones(Ng,1)*[2 1 2];
    end
end

TDP   = spm_MDP_VB_X(MDP,OPTIONS);
for i = 1:size(TDP,2)
    subplot(4,3,6 + i)
    spm_questions_plot(TDP(:,i))
end


return


function spm_questions_plot(MDP)
%% illustrate beliefs
%--------------------------------------------------------------------------

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
% label.modality{1} = 'noun';   label.outcome{1}  = {'square','triangle'};
% label.modality{2} = 'adj.';   label.outcome{2}  = {'green','red'};
% label.modality{3} = 'adverb'; label.outcome{3}  = {'above','below'};
% label.modality{4} = 'syntax'; label.outcome{4}  = {'1','2','3','Y','N','Sil.'};

% plot question, answer and posterior beliefs
%==========================================================================
cla;

for m = 1:numel(MDP)
    
    % Assemble question-and-answer
    %----------------------------------------------------------------------
    question = MDP(m).o(4,2);
    answer   = MDP(m).o(4,3);
    
    try
        qstr = cell2mat(MDP(m).MDP.label.outcome{1}(MDP(m).mdp(2).o));
        astr = cell2mat(MDP(m).MDP.label.outcome{1}(MDP(m).mdp(3).o));
    catch
        noun     = {'square','triangle'};
        adj      = {'green','red'};
        adverb   = {'above','below'};
        if question == 1
            qstr = ['is there a ' noun{MDP(m).o(1,2)} ' ?'];
        elseif question == 2
            qstr = ['Is any ' noun{MDP(m).o(1,2)} ' ' adverb{MDP(m).o(3,2)} ' ?'];
        elseif question == 3
            qstr = ['Is a ' adj{MDP(m).o(2,2)} ' ' noun{MDP(m).o(1,2)} ' ' adverb{MDP(m).o(3,2)} ' ?'];
        else
            qstr = '!';
        end
        if answer == 4
            astr = 'Yes !';
        elseif answer == 5
            astr = 'No !';
        else
            astr = 'I''m not sure';
        end
    end

    
    % is the answer right (for a single player)?
    %----------------------------------------------------------------------
    ind      = num2cell(MDP(m).s(:,3));
    if answer == find(MDP(m).A{4}(:,ind{:}),1) || m > 1
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
    col{1}  = [MDP(m).X{3}(2,T) MDP(m).X{3}(1,T) 0];
    col{1}  = col{1}*MDP(m).X{5}(1,T) + (1 - MDP(m).X{5}(1,T));
    col{2}  = [MDP(m).X{3}(2,T) MDP(m).X{3}(1,T) 0];
    col{2}  = col{2}*MDP(m).X{5}(2,T) + (1 - MDP(m).X{5}(2,T));
    col{3}  = [MDP(m).X{4}(2,T) MDP(m).X{4}(1,T) 0];
    col{3}  = col{3}*MDP(m).X{6}(1,T) + (1 - MDP(m).X{6}(1,T));
    col{4}  = [MDP(m).X{4}(2,T) MDP(m).X{4}(1,T) 0];
    col{4}  = col{4}*MDP(m).X{6}(2,T) + (1 - MDP(m).X{6}(2,T));
    
    plot(m,0,'^','MarkerSize',24,'LineWidth',4,'Color',col{4}), hold on
    plot(m,0,'s','MarkerSize',24,'LineWidth',4,'Color',col{3})
    plot(m,1,'^','MarkerSize',24,'LineWidth',4,'Color',col{2})
    plot(m,1,'s','MarkerSize',24,'LineWidth',4,'Color',col{1})
    
end

try, nq = MDP(1).n(4,2); catch, nq = 1; end
try, na = MDP(1).n(4,3); catch, na = 1; end
if ~nq, nq = 1; end
if ~na, na = 1; end
text(nq, 2,qstr,'HorizontalAlignment','Center')
text(na,-1,astr,'HorizontalAlignment','Center','FontWeight','bold','Color',cor)
axis([0 (m + 1) -1.5 2.5]), axis off, axis square


% return if multiple agents
%--------------------------------------------------------------------------
if numel(MDP) > 1, return, end

% upper and lower object
%--------------------------------------------------------------------------
rgb     = {[0 1 0],[1 0 0]};
shape   = {'s','^'};
plot(1 + 1/2,1,shape{MDP.s(5)},'MarkerSize',8,'LineWidth',1,'Color',rgb{MDP.s(3)})
plot(1 + 1/2,0,shape{MDP.s(6)},'MarkerSize',8,'LineWidth',1,'Color',rgb{MDP.s(4)})





 