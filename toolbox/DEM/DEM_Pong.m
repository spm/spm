function MDP = DEM_Pong
% Demo of planning as sophisticated inference (Tower of Hanoi)
%__________________________________________________________________________
%
% This routine demonstrates the acquisition of self organised behaviour
% using a simulated game of ping-pong. It starts by illustrating the
% emergence of game play during the acquisition of a likelihood model under
% precise beliefs about transitions through latent states. These
% transitions are generic and describe a (long) orbit with circular
% boundary conditions. This allows the agent to map successive states to
% observable outcomes. In this simulation the outcomes can be regarded as
% pixels (that can be on or off). In the particular set up here, a ball
% bounces around inside a box, while a paddle can be moved at the bottom
% right or left. Whenever there is a failure to hit the ball, unpredictable
% stimuli are delivered to all outcomes. This reflects an ambiguity that is
% specific to latent states under which the ball does not fall on the
% paddle. This ambiguity is, effectively, aversive; leading the agent to
% return the ball for short rallies, after it has acquired the requisite
% likelihood mapping (i.e., learned about the ambiguous outcomes).
% 
% We then repeat the simulation but making the box wider and therefore
% making the ball more difficult to track. This slows down learning and
% subsequent skilled gameplay. In the final simulation, this is rescued
% using inductive inference. With inductive constraints on active
% inference, skilled performance becomes fluent after about 512 timesteps
% (i.e., about two minutes of simulated time).
% 
% The depth of planning in the simulations is set to N = 1; i.e., looking
% ahead to time steps. The agents have mild preferences for not being
% overly stimulated; i.e., C{g} = [0;1].
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_dSprites.m 8447 2023-06-17 16:32:24Z karl $

%% set up and preliminaries
%==========================================================================
rng(1)

% size of problem
%--------------------------------------------------------------------------
Nr        = 6;                             % number of rows
Nc        = 5;                             % number of columns
[mdp,hid] = spm_Pong_MDP(Nr,Nc);

%% Simulate learning and subsequent performance
%==========================================================================

% Uninformative prior Dirichlet parameters
%--------------------------------------------------------------------------
[Nf,Ns,Nu,Ng] = spm_MDP_size(mdp);
for g = 1:Ng
    a{g} = spm_zeros(mdp.A{g}) + 1/16;
end

% specify a game and solve
%--------------------------------------------------------------------------
MDP   = mdp;
MDP.a = a;                                 % likelihood parameters
MDP.T = 512;                               % duration of game     
PDP   = spm_MDP_VB_XXX(MDP);

% illustrate Bayesian belief updating
%--------------------------------------------------------------------------
spm_figure('GetWin','Inference'); clf
spm_MDP_VB_trial(PDP);

% Illustrate Learning
%--------------------------------------------------------------------------
spm_figure('GetWin','Learning'); clf
spm_MDP_params(PDP)

% illustrate action and action selection
%--------------------------------------------------------------------------
spm_figure('GetWin','optimal behaviour'); clf
spm_report(PDP,Nr,Nc,hid)


%% Now repeat with a slightly more difficult problem
%==========================================================================

% size of problem
%--------------------------------------------------------------------------
Nr        = 4;                              % number of rows
Nc        = 8;                              % number of columns
[mdp,hid] = spm_Pong_MDP(Nr,Nc);

% Uninformative prior Dirichlet parameters
%--------------------------------------------------------------------------
[Nf,Ns,Nu,Ng] = spm_MDP_size(mdp);
for g = 1:Ng
    a{g} = spm_zeros(mdp.A{g}) + 1/16;
end

% specify a game and solve
%--------------------------------------------------------------------------
MDP   = mdp;
MDP.a = a;                                 % likelihood parameters
MDP.T = 512;                               % duration of game
PDP   = spm_MDP_VB_XXX(MDP);

% Illustrate Learning
%--------------------------------------------------------------------------
spm_figure('GetWin','Learning (2)'); clf
spm_MDP_params(PDP)

% illustrate action and action selection
%--------------------------------------------------------------------------
spm_figure('GetWin','Abductive inference'); clf
spm_report(PDP,Nr,Nc,hid)


%% Now repeat with inductive inference
%==========================================================================
MDP.hid = hid;
PDP     = spm_MDP_VB_XXX(MDP);

% illustrate action and action selection
%--------------------------------------------------------------------------
spm_figure('GetWin','Inductive inference'); clf
spm_report(PDP,Nr,Nc,hid)

% Update learned likelihood mapping and continue
%--------------------------------------------------------------------------
MDP   = spm_MDP_VB_update(MDP,PDP);
PDP   = spm_MDP_VB_XXX(MDP);

% illustrate action and action selection
%--------------------------------------------------------------------------
spm_figure('GetWin','Fluent behaviour'); clf
spm_report(PDP,Nr,Nc,hid)


return


% subroutines
%==========================================================================

function [MDP, hid] = spm_Pong_MDP(Nr,Nc)
% Creates an MDP structure for this game
%--------------------------------------------------------------------------
% Nr    = 6;                             % number of rows
% Nc    = 8;                             % number of columns
%
% hid   - Hidden states corresponding to hits
%--------------------------------------------------------------------------
Ng    = Nr*Nc;                           % number of locations
Ns    = 256;                             % maximum number of hidden states
for g = 1:Ng                             % default likelihood mapping
    A{g}        = zeros(2,Ns,Nc);
    A{g}(2,:,:) = 1;
end
B{1}  = zeros(Ns,Ns);                    % transition matrices

% Likelihood and transition tensors
%--------------------------------------------------------------------------
S     = zeros(Ng,4);

i     = 2;                               % Initial location (horizontal)
j     = 2;                               % Initial location (vertical)
p     = 1;                               % Momentum (horizontal)
q     = 1;                               % Momentum (vertical)
for s = 1:Ns

    % Check whether this state has been previously visited
    %----------------------------------------------------------------------
    r      = [i,j,p,q];
    k      = ismember(S,r,'rows');
    if any(k)
        r  = find(k);
        s  = s - 1;

        B{1}(s + 1,s,1) = 0;
        B{1}(r,s,1)     = 1;
        B{1}(r,s,1)     = 1;
        B{1} = B{1}(1:s,1:s);
        for g = 1:Ng
            A{g} = A{g}(:,1:s,:);
        end
        break
    else
        S(s,:) = r;
    end

    % Index of latent state
    %----------------------------------------------------------------------
    n  = sub2ind([Nr,Nc],i,j);

    A{n}(1,s,:)      = 1;
    A{n}(2,s,:)      = 0;
    B{1}(s + 1,s,1)  = 1;

    % uncomment to show orbit
    %----------------------------------------------------------------------
    % subplot(2,1,1)
    % plot(j,i,'o'), hold on
    % axis([1 Nc 1 Nr])
    % axis image, drawnow

    % Boundary conditions (switch momentum)
    %----------------------------------------------------------------------
    if ismember(i,[1,Nr]), p = -p; end
    if ismember(j,[1,Nc]), q = -q; end
    i  = i + p;
    j  = j + q;

end


% paddle (three actions)
%--------------------------------------------------------------------------
Nu    = 3;
B{2}  = zeros(Nc,Nc,Nu);
for u = 1:Nu
    B{2}(:,:,u) = spm_speye(Nc,Nc,u - 2,2);
end

% Add paddle to likelihood mapping (observations)
%--------------------------------------------------------------------------
for s = 1:Nc
    n = sub2ind([Nr,Nc],1,s);
    A{n}(1,:,s) = 1;
    A{n}(2,:,s) = 0;
end

% Enumerate the states and paths of the ensuing generative model
%--------------------------------------------------------------------------
Nf    = numel(B);                    % number of hidden factors
Ng    = numel(A);                    % number of outcome modalities
for f = 1:Nf
    Ns(f) = size(B{f},1);
    Nu(f) = size(B{f},3);
end
for g = 1:Ng
    No(g) = size(A{g},1);
end


% priors: (cost) C: mild preference for not being stimulated
%--------------------------------------------------------------------------
for g = 1:Ng
    C{g} = [0; 1];
end

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify prior beliefs about
% initial states (D) and paths through those states (E)
%--------------------------------------------------------------------------
for f = 1:Nf
    D{f} = sparse(1,1,1,Ns(f),1);     % First state
    E{f} = sparse(1,1,1,Nu(f),1);     % First path
    H{f} = ones(Ns(f),1);             % No intentions at this stage
end

% Ambiguity and latent states corresponding to hits
%--------------------------------------------------------------------------
hid    = [];
for s1 = 1:Ns(1)
    for s2 = 1:Ns(2)

        % Render misses ambiguous
        %------------------------------------------------------------------
        if S(s1,1) == 1 && s2 ~= S(s1,2)
            for g = 1:Ng
                A{g}(:,s1,s2) = 1;
            end
        end

        % And record latent states corresponding to hits
        %------------------------------------------------------------------
        if S(s1,1) == 1 && s2 == S(s1,2)
            hid(:,end + 1) = [s1;s2];
        end

    end
end

% specify controllable factors; here, the second factor
%--------------------------------------------------------------------------
U     = [0,1];                        % controllable factors

% Assemble MDP structure, with generative process
%==========================================================================
MDP.T = 8;                            % numer of moves
MDP.U = U;                            % controllable factors
MDP.A = A;                            % likelihood probabilities
MDP.B = B;                            % transition probabilities
MDP.C = C;                            % prior preferences
MDP.D = D;                            % prior over initial states
MDP.H = H;                            % prior over final states
MDP.E = E;                            % prior over initial paths
MDP.N = 1;                            % planning depth (2)

return


function spm_report(MDP,Nr,Nc,hid)
% Plots the inferred sequence of moves
%--------------------------------------------------------------------------
% If there are more than eight moves this subroutine will plot a movie

%% show sequence of moves
%==========================================================================
for t = 1:MDP.T

    % movie
    %----------------------------------------------------------------------
    if MDP.T > 8
        subplot(4,3,2), cla
    else
        subplot(4,4,t), cla
    end

    % plot
    %----------------------------------------------------------------------
    imagesc(1 - reshape(MDP.o(:,t),Nr,Nc)), axis image, axis xy
    title(sprintf('Time %i',t),'FontSize',12)
    drawnow

    % save movie
    %----------------------------------------------------------------------
    if MDP.T > 8
        I(t) = getframe(gca);
    end

end

% hits
%--------------------------------------------------------------------------
x  = find(ismember(MDP.s',hid','rows'));

subplot(6,1,3)
plot(MDP.F), hold on
plot(x,x*0 + mean(MDP.F),'.k','MarkerSize',16), hold off 
title('Negative variational free energy (ELBO)','FontSize',12)
xlabel('time'), ylabel('nats'), spm_axis tight

subplot(6,1,4)
plot(MDP.w), hold on
plot(x,x*0,'.k','MarkerSize',16), hold off 
title('Precision (confidence)','FontSize',12)
xlabel('time'), ylabel('nats'), spm_axis tight

subplot(6,1,5)
imagesc(MDP.R) 
title('Bayesian beliefs about policies','FontSize',12)
xlabel('time'), ylabel('policy')

% Place movie in graphic subject
%--------------------------------------------------------------------------
if MDP.T > 8
    subplot(4,3,2)
    set(gca,'Userdata',{I,32})
    set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')
end

return
