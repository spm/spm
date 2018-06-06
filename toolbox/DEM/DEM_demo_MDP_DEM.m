function MDP = DEM_demo_MDP_DEM
% Demo of mixed continuous and discrete state space modelling
%__________________________________________________________________________
%
% This routine  illustrates the combination of discrete and continuous
% state space models for active inference. In this example, the lowest
% level of a hierarchical Markov decision process (used to illustrate
% evidence accumulation during reading in related simulations) is equipped
% with a continuous time and state space dynamical model at the lowest
% level. This allows one to model both the categorical belief updates
% using belief propagation and the continuous belief updates using
% Bayesian filtering within the same model and associated inversion
% scheme.
%
% The key contribution of this  scheme is the message passing or belief
% propagation between the lowest discrete state (MDP) level and the
% highest level of the continuous state (DCM) models. In brief, during
% inversion, posterior beliefs about hidden causes of observable
% (continuous) inputs provide (probabilistic or posterior) outcomes for the
% (categorical) MDP scheme. In return, the posterior predictive density
% over outcomes of the MDP scheme specify priors on the hidden causes. In
% this example, these priors determine the salient locations to which the
% synthetic agent saccades. These saccades sample discriminative visual
% information that resolves uncertainty about the content of the local
% visual scene. Posterior expectations about the content then play the role
% of observations for higher (categorical) levels.
%
% Note that the priors from the MDP levels are time invariant (i.e., the
% attracting location of the saccade does not change during each saccadic
% epoch). Similarly, the posterior beliefs are over attributes that do
% not change during the saccadic sampling (i.e., the hidden cause of
% visual input at the attracting location). This underwrites a separation
% of temporal scales that is recapitulated at higher levels of the
% categorical model. The implementation of these schemes is as general as
% we could make it. The code below illustrates how one links MDP schemes
% to DPM schemes in a generic way through  hidden causes.
%
% More details about each level of the model are provided in line as
% annotated descriptions.
%
% see also: spm_MPD_VB_X.m
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_MDP_DEM.m 7326 2018-06-06 12:16:40Z karl $

% set up and preliminaries: first level
%--------------------------------------------------------------------------
rng('default')
pth = fileparts(mfilename('fullpath'));


% generative model at the sensory level (DEM): continuous states
%==========================================================================
% This level of model specification concerns the sampling of continuous
% data; here, visual stimuli encoded in a global structure (STIM). This is
% a generic specification that allows one to place various stimuli in the
% visual field – with user-specified parameters for sampling. The hidden
% causes of this generative model correspond to one attracting location
% and the content or stimulus that will be sampled at that location. The
% dynamics or hidden states in this level of the model are simple: the
% attracting location simply attracts the point of foveal fixation.


% memory map images and get hypotheses
%--------------------------------------------------------------------------
global STIM
str   = {'null','bird','seed','cats'};
for i = 1:numel(str)
    STIM.H{i} = spm_vol(fullfile(pth,[str{i},'.nii']));
end

% set-up:
%--------------------------------------------------------------------------
dim    = 32;                                  % dimension of visual sample
nb     = 6;                                   % number of basis functions
ns     = nb*nb;                               % number of sensory channels
STIM.W = 1/2;                                 % foveal sampling width
STIM.R = ones(dim,dim);                       % retinal precision
STIM.B = spm_dctmtx(dim,nb);                  % Basis functions (CRFs)
STIM.A = 4;                                   % Attenuation (spotlight)

% and locations
%--------------------------------------------------------------------------
L{1}   = [-8;-8];
L{2}   = [ 8;-8];
L{3}   = [-8; 8];
L{4}   = [ 8; 8];
STIM.L = L;


% mapping from outputs of higher (discrete) level to (hidden) causes
%==========================================================================

% true causes (U) and priors (C) for every combination of discrete states
%--------------------------------------------------------------------------
N     = 24;                                  % length of data sequence
nh    = length(STIM.H);                      % number of hypotheses
nl    = length(STIM.L);                      % number of locations
for i = 1:nh
    for j = 1:nl
        c           = [STIM.L{j}; sparse(i,1,1,nh,1)];
        u           = [STIM.L{j}; sparse(i,1,1,nh,1)];
        demi.U{i,j} = u*ones(1,N);
        demi.C{i,j} = c*ones(1,N);
    end
end

% evaluate true and priors over causes given discrete states
%--------------------------------------------------------------------------
o     = [4,1];
O{1}  = spm_softmax(sparse(1:4,1,1,nh,1));
O{2}  = spm_softmax(sparse(1,1,4,nl,1));

% generative model
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 2;                                 % order of
M(1).E.d = 1;                                 % generalised motion

% hidden states
%--------------------------------------------------------------------------
x      = [0;0];                               % oculomotor angle
v.x    = [8;0];                               % fixed (attracting) point
v.h    = sparse(nh,1);                        % hypothesis
g      = @ADEM_sample_image;

% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f = @(x,v,P) (v.x - x)/8;
M(1).g = @(x,v,P) spm_vec(x,g(x - v.x,v.h));
M(1).x = x;                                   % hidden states
M(1).V = 8;                                   % error precision (g)
M(1).W = 8;                                   % error precision (f)

% level 2:
%--------------------------------------------------------------------------
M(2).v = v;                                   % priors
M(2).V = [exp(8) exp(8) ones(1,nh)];


% generative process
%==========================================================================

% first level
%--------------------------------------------------------------------------
G(1).f = @(x,v,a,P) a;
G(1).g = @(x,v,a,P) spm_vec(x,g(x - v.x,v.h));
G(1).x = x;                                  % hidden states
G(1).V = exp(16);                            % error precision
G(1).W = exp(16);                            % error precision
G(1).U = [1 1 zeros(1,ns)];                  % gain

% second level
%--------------------------------------------------------------------------
G(2).v = v;                                  % exogenous forces
G(2).a = [0;0];                              % action forces
G(2).V = exp(16);

% generate and invert
%==========================================================================
DEM.G  = G;
DEM.M  = M;

% solve and save saccade
%--------------------------------------------------------------------------
DEM    = spm_MDP_DEM(DEM,demi,O,o);

% arrays (of functions) for graphics (in spm_MDP_plot)
%--------------------------------------------------------------------------
DEM.label{1} = @(x,v,P) g(x - v.x,v.h);
DEM.label{2} = @(x,v,P) 1 - sparse(ceil(x(1)*2 + 16),ceil(x(2)*2 + 16),1,32,32);

% overlay true values on the results of Bayesian filtering
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU,DEM.pU)

% show movies of sampling of continuous observations at this level
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_dem_mdp_movie(DEM)
drawnow


%  second level (discrete: lexical)
%==========================================================================
% There are two outcome modalities (what and where), encoding one of 4 cues
% and one of 4 locations.  The hidden states have four factors;
% corresponding to context (three words), where (the 4 locations)
% and two further factors modelling invariance. There are no specific prior
% preferences at this level, which means behaviour is purely epistemic.
%--------------------------------------------------------------------------
label.factor     = {'what','where','flip','flip'};
label.name{1}    = {'flee','feed','wait'};
label.name{2}    = {'1','2','3','4'};
label.name{3}    = {'-','+'};
label.name{4}    = {'-','+'};

label.modality   = {'what','where'};
label.outcome{1} = {'null','bird','seed','cat'};
label.outcome{2} = {'1','2','3','4'};


% prior beliefs about initial states
%--------------------------------------------------------------------------
D{1} = [1 1 1]';           % what:     {'flee','feed','wait'}
D{2} = [1 0 0 0]';         % where:    {'1',...,'4'}
D{3} = [8 1]';             % flip(ud): {'no','yes'}
D{4} = [8 1]';             % flip(rl): {'no','yes'}

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
                
                % A{1} what: {'null','bird','seed','cat'}
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


% MDP structure for this level (and subordinate DEM level)
%--------------------------------------------------------------------------
mdp.T     = 3;                      % number of updates
mdp.A     = A;                      % observation model
mdp.B     = B;                      % transition probabilities
mdp.D     = D;                      % prior over initial states
mdp.DEM   = DEM;
mdp.demi  = demi;

mdp.label = label;
mdp.Aname = label.modality;
mdp.Bname = label.factor;
mdp.alpha = 64;
mdp.chi   = 1/32;

clear A B D

MDP = spm_MDP_check(mdp);
clear mdp


% Third level (discrete: semantic)
%==========================================================================
% There are three outcome modalities (what, where and feedback), encoding
% one of three cues (i.e., words) in one of 4 locations.  The hidden states
% have three factors; corresponding to context (one of six sentences),
% with a particular order of words and the four successive locations. The
% third hidden factor corresponds to the categorical decision. A priori,
% the agent prefers to solicit correct feedback.
%--------------------------------------------------------------------------
label.factor     = {'story','where','decision'};
label.name{1}    = {...
    'flee wait feed wait'
    'wait wait wait feed'
    'wait flee wait feed'
    'flee wait feed flee'
    'wait wait wait flee'
    'wait flee wait flee'};
label.name{2}    = {'1st','2nd','3rd','4th'};
label.name{3}    = {'null','happy','sad'};

label.modality   = {'picture','where','feedback'};
label.outcome{1} = {'flee','feed','wait'};
label.outcome{2} = {'1','2','3','4'};
label.outcome{3} = {'null','right','wrong'};


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

% control states B(2): where {'stay' or 'proceed'}
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
C{3}(2,:) =  0;                 % the agent expects to be right
C{3}(3,:) = -4;                 % and not wrong

% MDP structure for this level (and subordinate MDP level)
% including links from outcomes at the current level to states of level
% below. This  complete the specification of the mixed hierarchical model
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

mdp.label = label;
mdp.alpha = 64;
mdp       = spm_MDP_check(mdp);


% invert this model of sentence reading
%==========================================================================
MDP  = spm_MDP_VB_X(mdp);

% show belief updates (and behaviour) over trials
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_MDP_VB_trial(MDP);

% illustrate phase-precession and electrophysiological responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
spm_MDP_VB_LFP(MDP,[],1);
subplot(3,1,3), spm_MDP_search_plot(MDP)

% electrophysiological responses over MDP levels
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
spm_MDP_VB_ERP(MDP,1);
subplot(4,1,4), spm_MDP_search_plot(MDP)

% illustrate evidence accumulation and perceptual synthesis (movie)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 6'); clf
spm_MDP_search_percept(MDP)

% movie of expected states and outcomes at all levels
%--------------------------------------------------------------------------
subplot(2,1,2)
spm_MDP_plot(MDP)


return




function spm_MDP_search_plot(MDP)
% illustrate visual search graphically
%--------------------------------------------------------------------------

% locations on page and of page
%--------------------------------------------------------------------------
x = [0 0;0 1;1 0;1 1];
y = [1 0;2 0;3 0;4 0]*3;
r = [-1,1]/2;

% load images
%--------------------------------------------------------------------------
load MDP_search_graphics

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
i   = 1:(size(T,1) - 8);
plot(T(i,1),T(i,2),'b ','LineWidth',1)
plot(X(:,1),X(:,2),'r.','MarkerSize',16)


function spm_MDP_search_percept(MDP)
% illustrates visual expectations with a movie
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

% plot posterior beliefs
%--------------------------------------------------------------------------
M     = [];
for p = 1:numel(MDP.mdp)
    mdp   = MDP.mdp(p);
    for k = 1:size(mdp.xn{1},4)
        for i = 1:size(mdp.xn{1},1)
            
            % movie over peristimulus time: first level expectations
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
                delete(hl);
            end
            for j = 1:numel(S)
                hl(j) = imagesc(r + y(MDP.o(2,p),1) + x(j,1),r + y(MDP.o(2,p),2) + x(j,2),S{j}/max(S{j}(:)));
            end
            
            
            % saccade
            %--------------------------------------------------------------
            X = y(MDP.o(2,p),:) + x(MDP.mdp(p).o(2,k),:);
            plot(X(:,1),X(:,2),'r.','MarkerSize',32)
            
            
            % movie over peristimulus time: second level expectations
            %--------------------------------------------------------------
            for j = 1:4
                S{j} = zeros(size(bird));
            end
            for f1 = 1:6
                
                % sequence of pictures for each story
                %----------------------------------------------------------
                if f1 == 1, a = {'cats','null','seed','null'}; end  % happy
                if f1 == 2, a = {'null','null','null','seed'}; end  % happy
                if f1 == 3, a = {'null','cats','null','seed'}; end  % happy
                if f1 == 4, a = {'cats','null','seed','cats'}; end  % sad
                if f1 == 5, a = {'null','null','null','cats'}; end  % sad
                if f1 == 6, a = {'null','cats','null','cats'}; end  % sad
                
                % mixture
                %----------------------------------------------------------
                for j = 1:numel(a)
                    S{j} = S{j} + eval(a{j})*MDP.xn{1}(end,f1,p,p);
                end
                
            end
            
            % image
            %--------------------------------------------------------------
            if i > 1
                delete(hh);
            end
            for j = 1:numel(S)
                hh(j) = imagesc(s + y(j,1) + 1/2,s + y(j,2) - 1 - 1/4,S{j}/max(S{j}(:)));
            end
            
            % categorisation
            %--------------------------------------------------------------
            if i > 1
                delete(ht);
            end
            pt    = MDP.xn{3}(end,2,p,p);
            ht(1) = text(12,-2.5,'happy','Color',[1,[1,1]*(1 - pt)],'FontSize',16);
            pt    = MDP.xn{3}(end,3,p,p);
            ht(2) = text(12,-2.0,'sad','Color',[[1,1]*(1 - pt),1],'FontSize',16);
            
            % save
            %--------------------------------------------------------------
            axis image ij, axis([2 14 -2 2]), drawnow
            if numel(M)
                M(end + 1) = getframe(gca);
            else
                M = getframe(gca);
            end
        end
    end
end

% save movie
%--------------------------------------------------------------------------
set(gca,'Userdata',{M,16})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
title('Narrative construction','FontSize',16)
