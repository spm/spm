function MDP = DEM_demo_MDP_reading
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
% how average behaviour (accuracy, number of saccades and saccade duration)
% depends upon prior preferences and prior precision.
%
% This demonstration uses a factorised version of the MDP scheme. In
% other words, we assume a mean field approximation to the posterior over
% different hidden states (context, location, scene reflection) – and over
% multiple modalities (what versus where).  This provides a parsimonious
% representation of posterior beliefs over hidden states – but does induce
% degree of overconfidence associated with approximate Bayesian inference.
%
% see also: DEM_demo_MDP_habits.m and spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_MDP_reading.m 6832 2016-07-11 09:28:44Z karl $

% set up and preliminaries
%==========================================================================
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
% rng('default')


% first level (lexcial)
%==========================================================================

% prior beliefs about initial states (in terms of counts_: D and d
%--------------------------------------------------------------------------
D{1} = [1 1 1]';           % what:     {'flee','feed','wait'}
D{2} = [1 0 0 0]';         % where:    {'1',...,'4'}
D{3} = [1 1]';             % flip(ud): {'no','yes'}
D{4} = [1 1]';             % flip(rl): {'no','yes'}


% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
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
                
                % A{1} what: {'null','bird,'seed','cat'}
                %==========================================================
                
                % saccade to cue location
                %----------------------------------------------------------
                A{1}(1,f1,f2,f3,f4) = strcmp(a{f2},'null');
                A{1}(2,f1,f2,f3,f4) = strcmp(a{f2},'bird');
                A{1}(3,f1,f2,f3,f4) = strcmp(a{f2},'seed');
                A{1}(4,f1,f2,f3,f4) = strcmp(a{f2},'cat');
                
                % A{2} where: {'1',...,'4'}
                %----------------------------------------------------------
                A{2}(f2,f1,f2,f3,f4) = 1;
                
            end
        end
    end
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


% MDP Structure
%--------------------------------------------------------------------------
mdp.T = 6;                      % number of moves
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.D = D;                      % prior over initial states

mdp.Aname = {'what','where'};
mdp.Bname = {'what','where','flip','flip'};

clear A B D

MDP = spm_MDP_check(mdp);
clear mdp

% second level (semantic)
%==========================================================================

% prior beliefs about initial states (in terms of counts_: D and d
%--------------------------------------------------------------------------
D{1} = [1 1 1 1 1 1]';   % what:   {'story 1',...,'story 6'}
D{2} = [1 0 0 0]';       % where:  {'1',...,'4'}
D{3} = [1 0 0]';         % report: {'null','happy','sad'}

% probabilistic mapping from hidden states to outcomes: A
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end
for f1 = 1:Ns(1)
    for f2 = 1:Ns(2)
        for f3 = 1:Ns(3)
            
            % sequence of pictures for each story
            %--------------------------------------------------------------
            if f1 == 1, a = {'flee','wait','feed','wait'}; end  % happy
            if f1 == 2, a = {'wait','wait','wait','feed'}; end  % happy
            if f1 == 3, a = {'wait','flee','wait','feed'}; end  % happy
            if f1 == 4, a = {'flee','wait','feed','flee'}; end  % sad
            if f1 == 5, a = {'wait','wait','wait','flee'}; end  % sad
            if f1 == 6, a = {'wait','flee','wait','flee'}; end  % sad
            
            
            % A{1} picture: 'flee','feed','wait'
            %==============================================================
            A{1}(1,f1,f2,f3) = strcmp(a{f2},'flee');
            A{1}(2,f1,f2,f3) = strcmp(a{f2},'feed');
            A{1}(3,f1,f2,f3) = strcmp(a{f2},'wait');
            
            % A{2} where: {'1',...,'4'}
            %--------------------------------------------------------------
            A{2}(f2,f1,f2,f3) = 1;
            
            % A{3} feedback: {'null','right','wrong'}
            %--------------------------------------------------------------
            hap = any(ismember([1 2 3],f1));
            sad = any(ismember([4 5 6],f1));
            A{3}(1,f1,f2,f3) = (f3 == 1);                         % undecided
            A{3}(2,f1,f2,f3) = (f3 == 2 & hap) | (f3 == 3 & sad); % right
            A{3}(3,f1,f2,f3) = (f3 == 3 & hap) | (f3 == 2 & sad); % wrong
            
        end
    end
end
Ng    = numel(A);
for g = 1:Ng
    No(g) = size(A{g},1);
end

% controlled transitions: B{f} for each factor
%--------------------------------------------------------------------------
for f = 1:Nf
    B{f} = eye(Ns(f));
end

% control states B(2): where {'stay','forward,'backward'}
%--------------------------------------------------------------------------
B{2}(:,:,1) = spm_speye(Ns(2),Ns(2), 0);
B{2}(:,:,2) = spm_speye(Ns(2),Ns(2),-1); B{2}(end,end,2) = 1;

% control states B(3): report {'null,'happy','sad'}
%--------------------------------------------------------------------------
for k = 1:Ns(3)
    B{3}(:,:,k) = 0;
    B{3}(k,:,k) = 1;
end

% allowable policies (specified as the next action) U
%--------------------------------------------------------------------------
U(1,1,:)  = [1 2 1]';           % move to next page
U(1,2,:)  = [1 1 2]';           % stay on current page and report happy
U(1,3,:)  = [1 1 3]';           % stay on current page and report sad

% priors: (utility) C
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g}  = zeros(No(g),1);
end
C{3}(2,:) =  2;                 % the agent expects to be right
C{3}(3,:) = -8;                 % and not wrong


% MDP Structure
%--------------------------------------------------------------------------
mdp.MDP  = MDP;
mdp.link = sparse(1,1,1,numel(MDP.D),Ng);

mdp.T = 5;                      % number of moves
mdp.U = U;                      % allowable policies
mdp.A = A;                      % observation model
mdp.B = B;                      % transition probabilities
mdp.C = C;                      % preferred outcomes
mdp.D = D;                      % prior over initial states
mdp.s = [1 1 1]';               % initial state

mdp.Aname = {'picture','where','feedback'};
mdp.Bname = {'story','where','decision'};

mdp  = spm_MDP_check(mdp);


% illustrate a single trial
%==========================================================================
MDP  = spm_MDP_VB_X(mdp);

% show belief updates (and behaviour)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_trial(MDP);

% illustrate phase-precession and responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],1);

subplot(3,1,3)
spm_MDP_search_plot(MDP)

% illustrate evidence accumulation and perceptual synthesis
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_search_percept(MDP)

return

% illustrate a sequence of trials
%==========================================================================

% true initial states – with context change at trial 12
%--------------------------------------------------------------------------

save

return


function spm_MDP_search_plot(MDP)
% illustrates visual search graphically
%--------------------------------------------------------------------------

% locations on page and of page
%--------------------------------------------------------------------------
x = [0 0;0 1;1 0;1 1];
y = [1 0;2 0;3 0;4 0]*3;
r = [-1,1]/2;

% load images
%--------------------------------------------------------------------------
load MDP_search_graphics
null = zeros(size(bird)) + 1;


% plot cues
%--------------------------------------------------------------------------
cla;
X     = [];
for p = 1:length(MDP.mdp)
    
    % latent cues for this hidden state
    %----------------------------------------------------------------------
    f    = MDP.mdp(p).s(:,1);
    
    if f(1) == 1, a = {'bird','cats';'null','null'}; end
    if f(1) == 2, a = {'bird','seed';'null','null'}; end
    if f(1) == 3, a = {'bird','null';'null','seed'}; end
    
    % flip cues according to hidden (invariants) states
    %----------------------------------------------------------------------
    if f(3) == 2, a = flipud(a); end
    if f(4) == 2, a = fliplr(a); end
    
    j     = MDP.s(2,p);
    for i = 1:numel(a)
        image(r + y(j,1) + x(i,1),r + y(j,2) + x(i,2), eval(a{i})), hold on
    end
    axis image ij, axis([2 14 -1 2])
    
    % Extract eye movements
    %----------------------------------------------------------------------
    for i = 1:numel(MDP.mdp(p).o(2,:))
        X(end + 1,:) = y(MDP.o(2,p),:) + x(MDP.mdp(p).o(2,i),:);
    end
end

% Smooth and plot eye movements
%--------------------------------------------------------------------------
for j = 1:2
    T(:,j) = interp(X(:,j),8,2);
    T(:,j) = T(:,j) + spm_conv(randn(size(T(:,j))),2)/16;
end
plot(T(:,1),T(:,2),'b ','LineWidth',1)
plot(X(:,1),X(:,2),'r.','MarkerSize',16)


function spm_MDP_search_percept(MDP)
% illustrates visual search graphically
%--------------------------------------------------------------------------
clf; 
subplot(4,1,1), spm_MDP_search_plot(MDP)
axis image ij, axis([2 14 -2 2]),
subplot(4,1,2)
axis image ij, axis([2 14 -2 2]), hold on

% locations on page and of page
%--------------------------------------------------------------------------
x = [0 0;0 1;1 0;1 1];
y = [1 0;2 0;3 0;4 0]*3;
r = [-1,1]/2;
s = r*(1 + 1/4);


% load images
%--------------------------------------------------------------------------
load MDP_search_graphics
null = zeros(size(bird)) + 1;

mdp = MDP.mdp(1);
try
    d = mdp.D;
catch
    d = mdp.d;
end
Nf    = numel(d);
for f = 1:Nf
    Ns(f) = numel(d{f});
end

% plot pposterior beliefs
%--------------------------------------------------------------------------
M     = [];
for p = 1:numel(MDP.mdp)
    mdp   = MDP.mdp(p);
    for k = 1:size(mdp.xn{1},4)
        for i = 1:size(mdp.xn{1},1)
            
            % movie over peristimulus time
            %--------------------------------------------------------------
            for j = 1:4
                S{j} = zeros(size(bird));
            end
            for f1 = 1:Ns(1)
                for f2 = 1:Ns(2)
                    for f3 = 1:Ns(3)
                        for f4 = 1:Ns(4)
                            
                            % latent cues for this hidden state
                            %----------------------------------------------
                            if f1 == 1, a = {'bird','cats';'null','null'}; end
                            if f1 == 2, a = {'bird','seed';'null','null'}; end
                            if f1 == 3, a = {'bird','null';'null','seed'}; end
                            
                            % flip cues according to hidden (invariants) states
                            %----------------------------------------------
                            if f3 == 2, a = flipud(a); end
                            if f4 == 2, a = fliplr(a); end
                            
                            % mixture
                            %----------------------------------------------
                            q     = mdp.xn{1}(i,f1,1,k)*mdp.xn{3}(i,f3,1,k)*mdp.xn{4}(i,f4,1,k);
                            for j = 1:4
                                S{j} = S{j} + eval(a{j})*q;
                            end
                        end
                    end
                end
            end
            
            % image
            %--------------------------------------------------------------
            if i > 1
                delete(h);
            end
            for j = 1:numel(S)
                h(j) = imagesc(r + y(MDP.o(2,p),1) + x(j,1),r + y(MDP.o(2,p),2) + x(j,2),S{j}/max(S{j}(:)));
            end
            
            
            % saccade
            %--------------------------------------------------------------
            X = y(MDP.o(2,p),:) + x(MDP.mdp(p).o(2,k),:);
            plot(X(:,1),X(:,2),'r.','MarkerSize',32)
            
            % save
            %--------------------------------------------------------------
            axis image ij, axis([2 14 -2 2]), drawnow
            if numel(M)
              %  M(end + 1) = getframe(gca);
            else
                M = getframe(gca);
            end
            
        end
        
    end
    for i = 1:size(MDP.xn{1},1)
        
        % movie over peristimulus time
        %------------------------------------------------------------------
        for j = 1:4
            S{j} = zeros(size(bird));
        end
        for f1 = 1:6
            
            % sequence of pictures for each story
            %--------------------------------------------------------------
            if f1 == 1, a = {'cats','null','seed','null'}; end  % happy
            if f1 == 2, a = {'null','null','null','seed'}; end  % happy
            if f1 == 3, a = {'null','cats','null','seed'}; end  % happy
            if f1 == 4, a = {'cats','null','seed','cats'}; end  % sad
            if f1 == 5, a = {'null','null','null','cats'}; end  % sad
            if f1 == 6, a = {'null','cats','null','cats'}; end  % sad
            
            % mixture
            %--------------------------------------------------
            for j = 1:numel(a)
                S{j} = S{j} + eval(a{j})*MDP.xn{1}(i,f1,1,p);
            end
            
        end
        
        % image
        %------------------------------------------------------------------
        if i > 1
            delete(h);
        end
        for j = 1:numel(S)
            h(j) = imagesc(s + y(j,1) + 1/2,s + y(j,2) - 1 - 1/4,S{j}/max(S{j}(:)));
        end
        
        
        % save
        %------------------------------------------------------------------
        drawnow
        % M(end + 1) = getframe(gca);
        
        
    end
end

% save movie
%--------------------------------------------------------------------------
set(gca,'Userdata',{M,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Narrative construction','FontSize',16)



