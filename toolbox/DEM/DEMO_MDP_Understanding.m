function mdp = DEMO_MDP_Understanding
% This demo uses a hierarchical version of the ubiquitous T-maze demo, in
% which the agent is equipped with a space of hypotheses about why it chose
% to act in a certain way. This means that, when queried, it is able to
% communicate an explanation for its actions.
%
% see also: DEM_demo_MDP_X.m and spm_MPD_VB_X.m
%__________________________________________________________________________

% Thomas Parr
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% set up and preliminaries
%==========================================================================
rng('default')

% FIRST LEVEL MODEL (T-MAZE)
%--------------------------------------------------------------------------

% Linguistic hidden states:
syn      = {'...', 'Please', 'explain', 'your', '#1', 'move', ' ', '-', 'I', '#2', 'the reward', 'was on the', '#3', 'so I', '#4', 'by going to', '#5'};
sem{1}   = {'first','second'};
sem{2}   = {'knew','did not know'};
sem{3}   = {'left','right'};
sem{4}   = {'guessed','explored','found it'};
sem{5}   = {'the left.','the right.','the cue.'};
words{1} = {' ', 'Please', 'explain', 'your', 'first', 'second', 'move'};
words{2} = {' ', 'I', 'knew', 'did not know', 'the reward', 'was on the', 'right', 'left' 'so I', 'guessed', 'explored', 'found it', 'by going to', 'the right.', 'the left.', 'the cue.'};

% outcome probabilities: A
%--------------------------------------------------------------------------
% We start by specifying the probabilistic mapping from hidden states
% to outcomes; where outcome can be exteroceptive or interoceptive: The
% exteroceptive outcomes A{1} provide cues about location and context,
% while interoceptive outcome A{2) denotes different levels of reward
%--------------------------------------------------------------------------
a      = .98;
b      = 1 - a;
for f3 = 1:17
    for f4 = 1:2
        for f5 = 1:2
            for f6 = 1:3
                for f7 = 1:3
                    if f3 == 1 % If solving the maze
                        A{1}(:,:,1,f3,f4,f5,f6,f7) = [...
                            1 0 0 0;    % cue start
                            0 1 0 0;    % cue left
                            0 0 1 0;    % cue right
                            0 0 0 1     % cue CS right
                            0 0 0 0];   % cue CS left
                        A{1}(:,:,2,f3,f4,f5,f6,f7) = [...
                            1 0 0 0;    % cue start
                            0 1 0 0;    % cue left
                            0 0 1 0;    % cue right
                            0 0 0 0     % cue CS right
                            0 0 0 1];   % cue CS left
                        
                        A{2}(:,:,1,f3,f4,f5,f6,f7) = [...
                            1 0 0 1;    % reward neutral
                            0 a b 0;    % reward positive
                            0 b a 0];   % reward negative
                        A{2}(:,:,2,f3,f4,f5,f6, f7) = [...
                            1 0 0 1;    % reward neutral
                            0 b a 0;    % reward positive
                            0 a b 0];   % reward negative
                        
                        % Audition (questioner)
                        A{3}(1,:,:,f3,f4,f5,f6, f7) = ones(1,4,2); % Null
                        
                        % Vocalisation
                        A{4}(1,:,:,f3,f4,f5,f6, f7) = ones(1,4,2); % Null
                        
                    else % If explaining the solution
                        A{1}(:,:,:,f3,f4,f5,f6,f7) = repmat([1;zeros(4,1)],1,4,2);
                        A{2}(:,:,:,f3,f4,f5,f6,f7) = repmat([1;zeros(2,1)],1,4,2);
                        
                        if f3 < 8
                            % Vocalisation
                            A{4}(1,:,:,f3,f4,f5,f6, f7) = ones(1,4,2); % Null
                            
                            for g = 1:length(words{1})
                                if syn{f3}(1) ~= '#'
                                    A{3}(g,:,:,f3,f4,f5,f6,f7) = strcmp(words{1}{g},syn{f3})*ones(1,4,2);
                                elseif syn{f3}(2) == '1'
                                    A{3}(g,:,:,f3,f4,f5,f6,f7) = strcmp(words{1}{g},sem{1}{f4})*ones(1,4,2);
                                end
                            end
                            
                        else
                            % Audition (questioner)
                            A{3}(1,:,:,f3,f4,f5,f6,f7) = ones(1,4,2); % Null
                            for f2 = 1:2
                                for g = 1:length(words{2})
                                    if syn{f3}(1) ~= '#'
                                        if syn{f3}(1) == '-'
                                            A{4}(g,:,f2,f3,f4,f5,f6,f7) = strcmp(words{2}{g},' ')*ones(1,4);
                                        else
                                            A{4}(g,:,f2,f3,f4,f5,f6,f7) = strcmp(words{2}{g},syn{f3})*ones(1,4);
                                        end
                                    elseif syn{f3}(2) == '2'
                                        A{4}(g,:,f2,f3,f4,f5,f6,f7) = strcmp(words{2}{g},sem{2}{f5})*ones(1,4);
                                    elseif syn{f3}(2) == '3'
                                        A{4}(g,:,f2,f3,f4,f5,f6,f7) = strcmp(words{2}{g},sem{3}{f2})*ones(1,4);
                                    elseif syn{f3}(2) == '4'
                                        A{4}(g,:,f2,f3,f4,f5,f6,f7) = strcmp(words{2}{g},sem{4}{f6})*ones(1,4);     
                                    elseif syn{f3}(2) == '5'
                                        A{4}(g,:,f2,f3,f4,f5,f6,f7) = strcmp(words{2}{g},sem{5}{f7})*ones(1,4);         
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

% controlled transitions: B{u}
%--------------------------------------------------------------------------
% Next, we have to specify the probabilistic transitions of hidden states
% for each factor. Here, there are four actions taking the agent directly
% to each of the four locations.
%--------------------------------------------------------------------------
B{1}(:,:,1)  = [1 0 0 1; 0 1 0 0;0 0 1 0;0 0 0 0];
B{1}(:,:,2)  = [0 0 0 0; 1 1 0 1;0 0 1 0;0 0 0 0];
B{1}(:,:,3)  = [0 0 0 0; 0 1 0 0;1 0 1 1;0 0 0 0];
B{1}(:,:,4)  = [0 0 0 0; 0 1 0 0;0 0 1 0;1 0 0 1];

% context/semantics, which cannot be changed by action
%--------------------------------------------------------------------------
B{2}  = eye(2);
B{4}  = eye(2);
B{5}  = eye(2);
B{6}  = eye(3);
B{7}  = eye(3);

% Syntactic states
%--------------------------------------------------------------------------
B{3}  = zeros(17);
B{3}(1,1)  = 1;
B{3}(3:7,2:6) = eye(5);
B{3}(7,7)  = 1;
B{3}(9:end,8:end-1) = eye(9);
B{3}(end,end) = 1;

% priors: (utility) C
%--------------------------------------------------------------------------
% Finally, we have to specify the prior preferences in terms of log
% probabilities over outcomes. Here, the agent prefers rewards to losses -
% and does not like to be exposed
%--------------------------------------------------------------------------
C{1}  = [repmat([-1;0;0;0;0],1,3),zeros(5,7)];
c     = 6;
C{2}  = [repmat([0;c;-c],1,3),zeros(3,7)];
C{3}  = zeros(size(A{3},1),10);
C{4}  = zeros(size(A{4},1),10);

% now specify prior beliefs about initial states, in terms of counts. Here
% the hidden states are factorised into location and context:
%--------------------------------------------------------------------------

%   Maze states
%--------------------------------------------------------------------------
D{1} = [1 0 0 0]';      % Location
D{2} = [1 1]'/2;        % Context (left or right)

%   Linguistic states
%--------------------------------------------------------------------------
D{3} = ones(17,1)/17;   % Syntactic states
D{4} = [1 1]'/2;        % #1 (first/second)
D{5} = [1 1]'/2;        % #2 (knew/didn't know)
D{6} = [1 1 1]'/3;      % #4 (guessed/explored/found it)
D{7} = [1 1 1]'/3;      % #5 (the left/the right/the cue)

% allowable policies (of depth T).  These are just sequences of actions
% (with an action for each hidden factor)
%--------------------------------------------------------------------------
V(:,:,1) = [1  1  1  1  2  3  4  4  4  4
    repmat([1  2  3  4  2  3  1  2  3  4],8,1)];
V(:,:,2) = 1;
V(:,:,3) = 1;
V(:,:,4) = 1;
V(:,:,5) = 1;
V(:,:,6) = 1;
V(:,:,7) = 1;

% Self-generated outcomes
%--------------------------------------------------------------------------
n = zeros(numel(A),size(V,1)+1);
n(4,:) = 1;

% MDP Structure - this will be used to generate arrays for multiple trials
%==========================================================================
mdp.V = V;                       % allowable policies
mdp.A = A;                       % observation model
mdp.B = B;                       % transition probabilities
mdp.C = C;                       % preferred outcomes
mdp.D = D;                       % prior over initial states
mdp.E = ones(10,1)/10;           % prior over policies
mdp.chi = -exp(64);              % preclude early termination
mdp.n = n;                       % outcome modality x time (1s for self generated outcomes, and 0s otherwise)
mdp.tau   = 12;                  % time constant for belief updating

% Solve lower level model (maze) to check behaviour:
%--------------------------------------------------------------------------
% MDP      = mdp;
% MDP.s    = ones(7,1);
% MDP.s(3) = 2;
% % MDP.D{3} = zeros(16,1); MDP.D{3}(8) = 1;
% MDP   = spm_MDP_VB_X(MDP);
% spm_figure('GetWin','Figure 1'); clf
% MDP_understanding_animation(MDP,words)

clear MDP

% SECOND LEVEL MODEL (EXPLAINING BEHAVIOUR)
%--------------------------------------------------------------------------
clear A B C D

% Higher level hidden state priors
%--------------------------------------------------------------------------
D{1} = [1 0 0]'/3;   % Narrative structure (solve maze, query, respond)
D{2} = [1 1 1 1]'/4; % First move
D{3} = [1 1 1 1]'/4; % Second move
D{4} = [1 1]'/2;     % Context
D{5} = [1 1]'/2;     % Query (first/second move)

% Higher level likelihood
%--------------------------------------------------------------------------
for f1 = 1:length(D{1})
    for f2 = 1:length(D{2})
        for f3 = 1:length(D{3})
            for f4 = 1:length(D{4})
                for f5 = 1:length(D{5})
% Initial syntactic state
%--------------------------------------------------------------------------
                    A{1}(1,f1,f2,f3,f4,f5) = f1 == 1;
                    A{1}(2,f1,f2,f3,f4,f5) = f1 == 2;
                    A{1}(8,f1,f2,f3,f4,f5) = f1 == 3;
                    A{1}(numel(syn),f1,f2,f3,f4,f5) = 0;

% Semantic state #1
%--------------------------------------------------------------------------
                    A{2}(f5,f1,f2,f3,f4,f5) = 1;
                    
% Semantic state #2
%--------------------------------------------------------------------------
                    if f5 == 1 % if first move is queried...
                        if f2 == 1 % ...and if stays put...
                            A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                        elseif f2 == 2 %...or if goes to left arm...
                            if f4 == 1 %...and this is the correct context...
                                A{3}(1,f1,f2,f3,f4,f5) = 1; %...conclude did know
                            else %...but if incorrect context...
                                A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                            end
                        elseif f2 == 3 %...or if goes to right arm...
                            if f4 == 2 %...and this is the correct context...
                                A{3}(1,f1,f2,f3,f4,f5) = 1; %...conclude did know
                            else %...but if incorrect context...
                                A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                            end
                        else %...or if goes to cue location...
                            A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                        end
                    else % if second move is queried...
                        if f3 == 1 % ...and if stays put...
                            A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                        elseif f3 == 2 %...or if goes to left arm...
                            if f4 == 1 %...and this is the correct context...
                                A{3}(1,f1,f2,f3,f4,f5) = 1; %...conclude did know
                            else %...but if incorrect context...
                                A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                            end
                        elseif f3 == 3 %...or if goes to right arm...
                            if f4 == 2 %...and this is the correct context...
                                A{3}(1,f1,f2,f3,f4,f5) = 1; %...conclude did know
                            else %...but if incorrect context...
                                A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                            end
                        else %...or if goes to cue location...
                            A{3}(2,f1,f2,f3,f4,f5) = 1; %...conclude did not know
                        end
                    end
% Semantic state #3
%--------------------------------------------------------------------------
                    A{4}(f4,f1,f2,f3,f4,f5) = 1;
                    
% Semantic state #4
%--------------------------------------------------------------------------
                    if f5 == 1 % if first move is queried...
                        if f2 == 1 % ...and if stays put...
                            A{5}(:,f1,f2,f3,f4,f5) = 1; %...no useful information
                        elseif f2 == 2 %...or if goes to left arm...
                            if f4 == 1 %...and this is the correct context...
                                A{5}(3,f1,f2,f3,f4,f5) = 1; %...found it
                            else %...but if incorrect context...
                                A{5}(1,f1,f2,f3,f4,f5) = 1; %...conclude guessed
                            end
                        elseif f2 == 3 %...or if goes to right arm...
                            if f4 == 2 %...and this is the correct context...
                                A{5}(3,f1,f2,f3,f4,f5) = 1; %...found it
                            else %...but if incorrect context...
                                A{5}(1,f1,f2,f3,f4,f5) = 1; %...conclude guessed
                            end
                        else %...or if goes to cue location...
                            A{5}(2,f1,f2,f3,f4,f5) = 1; %...explored
                        end
                    else % if second move is queried...
                        if f3 == 1 % ...and if stays put...
                            A{5}(:,f1,f2,f3,f4,f5) = 1; %...no useful information
                        elseif f3 == 2 %...or if goes to left arm...
                            if f4 == 1 %...and this is the correct context...
                                A{5}(3,f1,f2,f3,f4,f5) = 1; %...found it
                            else %...but if incorrect context...
                                A{5}(1,f1,f2,f3,f4,f5) = 1; %...conclude guessed
                            end
                        elseif f3 == 3 %...or if goes to right arm...
                            if f4 == 2 %...and this is the correct context...
                                A{5}(3,f1,f2,f3,f4,f5) = 1; %...found it
                            else %...but if incorrect context...
                                A{5}(1,f1,f2,f3,f4,f5) = 1; %...conclude guessed
                            end
                        else %...or if goes to cue location...
                            A{5}(2,f1,f2,f3,f4,f5) = 1; %...explored
                        end
                    end
% Semantic state #5
%--------------------------------------------------------------------------
                    if f5 == 1 % if first move is queried...
                        if f2 == 1 % ...and if stays put...
                            A{6}(:,f1,f2,f3,f4,f5) = 1; %...no useful information
                        elseif f2 == 2 %...or if goes to left arm...
                            A{6}(1,f1,f2,f3,f4,f5) = 1; %...left arm
                        elseif f2 == 3 %...or if goes to right arm...
                            A{6}(2,f1,f2,f3,f4,f5) = 1; %...right arm
                        else %...or if goes to cue location...
                            A{6}(3,f1,f2,f3,f4,f5) = 1; %...cue arm
                        end
                    else % if second move is queried...
                        if f2 == 2 % if first move was to left arm (absorbing state)
                            A{6}(1,f1,f2,f3,f4,f5) = 1; %...left arm
                        elseif f2 == 3 % if first move was to right arm (absorbing state)
                            A{6}(2,f1,f2,f3,f4,f5) = 1; %...right arm
                        elseif f3 == 1 % centre of maze
                            A{6}(:,f1,f2,f3,f4,f5) = 1; %...no useful information
                        elseif f3 == 2 %...or if goes to left arm...
                            A{6}(1,f1,f2,f3,f4,f5) = 1; %...left arm
                        elseif f3 == 3 %...or if goes to right arm...
                            A{6}(2,f1,f2,f3,f4,f5) = 1; %...right arm
                        else %...or if goes to cue location...
                            A{6}(3,f1,f2,f3,f4,f5) = 1; %...cue arm
                        end
                    end
% Policy
%--------------------------------------------------------------------------
                    if f2 == 2 || f2 == 3
                        A{7}(mdp.V(1,:,1)==f2,f1,f2,f3,f4,f5) = 1;
                    else
                        ind1 = find(mdp.V(1,:,1)==f2);
                        ind2 = find(mdp.V(2,:,1)==f3);
                        ind  = intersect(ind1,ind2);
                        A{7}(ind,f1,f2,f3,f4,f5) = 1;
                    end
                end
            end
        end
    end
end

% Transition probabilities
%--------------------------------------------------------------------------
B{1} = [0 0 0; 1 0 0; 0 1 1];
B{2} = eye(4);
B{3} = eye(4);
B{4} = eye(2);
B{5} = eye(2);

% Preferences
%--------------------------------------------------------------------------
for g = 1:numel(A)
    C{g} = ones(size(A{g},1),1);
end

% Construct MDP
%--------------------------------------------------------------------------
MDP.MDP = mdp; clear mdp

MDP.A = A;
MDP.B = B;
MDP.C = C;
MDP.D = D;
MDP.T = 3;
MDP.link       = zeros(numel(MDP.MDP.D),numel(MDP.A));
MDP.link(3,1)  = 1;
MDP.link(4,2)  = 1;
MDP.link(5,3)  = 1;
MDP.link(2,4)  = 1;
MDP.link(6,5)  = 1;
MDP.link(7,6)  = 1;
MDP.linkE      = zeros(1,numel(MDP.A));
MDP.linkE(end) = 1;

% For first move
%--------------------------------------------------------------------------
% MDP.s    = [1;1;1;2;1];  

% For confabulation
%--------------------------------------------------------------------------
% MDP.D{1} = [0 1 0]'/3;   
% MDP.s    = [2;1;1;2;1];  

mdp = spm_MDP_VB_X(MDP);

% Plot beliefs and behaviour (animation)
%--------------------------------------------------------------------------
spm_figure('GetWin','Animation'); clf
MDP_understanding_animation(mdp,words)

spm_figure('GetWin','Solve maze'); clf
MDP_understanding_maze(mdp.mdp(1))

spm_figure('GetWin','Query'); clf
MDP_understanding_maze(mdp.mdp(2))

spm_figure('GetWin','Answer'); clf
MDP_understanding_maze(mdp.mdp(3))

% Plot beliefs over time
%--------------------------------------------------------------------------
spm_figure('GetWin','Beliefs'); clf
MDP_understanding_beliefs(mdp)


function MDP_understanding_animation(MDP,words)
% If the input is a single level (inverted) MDP, this function plots the
% behaviour elicited from and inferences drawn by this model as an
% animation. If it is a hierarchical MDP, the function is applied
% successively to each of the lower level MDPs.

if isfield(MDP,'mdp')
    for i = 1:length(MDP.mdp)
        MDP_understanding_animation(MDP.mdp(i),words)
    end
else
    
% Maze structure
%--------------------------------------------------------------------------
    M = [0 0 0 0 0;
        0 1 1 1 0;
        0 0 1 0 0;
        0 0 1 0 0;
        0 0 0 0 0];
    
% Coordinates for each location
%--------------------------------------------------------------------------
    L = [3 3; 2 2; 4 2; 3 4 ;3 4];
    
    % Trajectory (with added noise)
    %--------------------------------------------------------------------------
    for t = 1:MDP.T
        X(t,:) = L(MDP.o(1,t),:);
    end
    for j = 1:2
        x(:,j) = interp(X(:,j),8,2);
    end
    x = x + spm_conv(randn(size(x)),2)/16;
    
% Animate
%--------------------------------------------------------------------------
    for t = 1:MDP.T
% Plot conversation
%--------------------------------------------------------------------------
        subplot(3,1,2)
        cla
        txt = [];
        for tt = 1:t
            txt = [txt ' ' words{1}{MDP.o(3,tt)}];
        end
        axis([0 10 -4 2])
        text(0,1,txt,'FontSize',12)
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        box on
        title('Query')
        hold off
        
        subplot(3,1,3)
        cla
        txt =  [];
        txt2 = [];
        for tt = 1:t
            if tt < 8
                txt = [txt ' ' words{2}{MDP.o(4,tt)}];
            else
                txt2 = [txt2 ' ' words{2}{MDP.o(4,tt)}];
            end
        end
        axis([0 10 -4 2])
        text(0,1,txt,'FontSize',12)
        if ~isempty(txt2)
        text(0,0,txt2,'FontSize',12)
        end
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        box on
        title('Response')
        hold off
        
        for i = 1:size(MDP.xn,1)
            
% Plot maze and associated beliefs
%--------------------------------------------------------------------------
            subplot(3,1,1)
            imagesc(M), colormap gray, axis equal, axis off, hold on
            title('Maze')
            
            % Plot beliefs about location
            for j = 1:4
                plot(L(j,1),L(j,2),'o','MarkerSize',16,'Color',[1,1 - MDP.xn{1}(i,j,t,t),1 - MDP.xn{1}(i,j,t,t)])
            end
            
            % Plot beliefs about context
            plot(L(2,1),L(2,2),'.','Markersize',30,'Color',[1 - MDP.xn{2}(i,1,t,t),1,1 - MDP.xn{2}(i,1,t,t)])
            plot(L(3,1),L(3,2),'.','Markersize',30,'Color',[1 - MDP.xn{2}(i,2,t,t),1,1 - MDP.xn{2}(i,2,t,t)])
            text(L(4,1),L(4,2),'\leftarrow','Color',[1 - MDP.xn{2}(i,1,t,t),1,1 - MDP.xn{2}(i,1,t,t)],'FontWeight','bold', 'HorizontalAlignment', 'center','Interpreter','Tex')
            text(L(4,1),L(4,2),'\rightarrow','Color',[1 - MDP.xn{2}(i,2,t,t),1,1 - MDP.xn{2}(i,2,t,t)],'FontWeight','bold', 'HorizontalAlignment', 'center','Interpreter','Tex')
            
            % Plot location
            plot(L(MDP.o(1,t),1),L(MDP.o(1,t),2),'.r','MarkerSize',16), hold off
            
% Plot movement
%--------------------------------------------------------------------------
            if t < MDP.T
                for k = (1+(t-1)*8):(1+t*8)
                    % Plot maze
                    subplot(3,1,1)
                    imagesc(M), colormap gray, axis equal, axis off, hold on
                    title('Maze')
                    
                    % Plot beliefs about location
                    for j = 1:4
                        plot(L(j,1),L(j,2),'o','MarkerSize',16,'Color',[1,1 - MDP.xn{1}(i,j,t,t),1 - MDP.xn{1}(i,j,t,t)])
                    end
                    
                    % Plot beliefs about context
                    plot(L(2,1),L(2,2),'.','Markersize',30,'Color',[1 - MDP.xn{2}(i,1,t,t),1,1 - MDP.xn{2}(i,1,t,t)])
                    plot(L(3,1),L(3,2),'.','Markersize',30,'Color',[1 - MDP.xn{2}(i,2,t,t),1,1 - MDP.xn{2}(i,2,t,t)])
                    text(L(4,1),L(4,2),'\leftarrow','Color',[1 - MDP.xn{2}(i,1,t,t),1,1 - MDP.xn{2}(i,1,t,t)],'FontWeight','bold', 'HorizontalAlignment', 'center','Interpreter','Tex')
                    text(L(4,1),L(4,2),'\rightarrow','Color',[1 - MDP.xn{2}(i,2,t,t),1,1 - MDP.xn{2}(i,2,t,t)],'FontWeight','bold', 'HorizontalAlignment', 'center','Interpreter','Tex')
                    
                    % Plot location
                    plot(x(k,1),x(k,2),'.r','MarkerSize',16), pause(0.1)
                    hold off
                end
            end
        end
    end
end
    

function MDP_understanding_maze(MDP)
% Maze structure
%--------------------------------------------------------------------------
M = [0 0 0 0 0;
    0 1 1 1 0;
    0 0 1 0 0;
    0 0 1 0 0;
    0 0 0 0 0];

% Coordinates for each location
%--------------------------------------------------------------------------
L = [3 3; 2 2; 4 2; 3 4 ;3 4];
for t = 1:3
    subplot(3,1,t)
    imagesc(M), colormap gray, axis equal, axis off, hold on
    title(['t = ' num2str(t)])
    
    % Plot beliefs about location
    xn  = MDP.xn{1}(end,:,t,t);
    xn(xn<0) = 0;
    xn(xn>1) = 1;
    for j = 1:4
        plot(L(j,1),L(j,2),'o','MarkerSize',16,'Color',[1,1 - xn(j),1 - xn(j)])
    end
    
    % Plot beliefs about context
    xn  = MDP.xn{2}(end,:,t,t);
    xn(xn<0) = 0;
    xn(xn>1) = 1;
    plot(L(2,1),L(2,2),'.','Markersize',30,'Color',[1 - xn(1),1,1 - xn(1)])
    plot(L(3,1),L(3,2),'.','Markersize',30,'Color',[1 - xn(2),1,1 - xn(2)])
    text(L(4,1),L(4,2),'\leftarrow','Color',[1 - xn(1),1,1 - xn(1)],'FontWeight','bold', 'HorizontalAlignment', 'center','Interpreter','Tex')
    text(L(4,1),L(4,2),'\rightarrow','Color',[1 - xn(2),1,1 - xn(2)],'FontWeight','bold', 'HorizontalAlignment', 'center','Interpreter','Tex')
   
    % Plot location
    plot(L(MDP.o(1,t),1),L(MDP.o(1,t),2),'.r','MarkerSize',16), hold off
end

function MDP_understanding_beliefs(MDP)
XX = [];
for i = 1:numel(MDP.xn)
    X = [];
    for j = 1:size(MDP.xn{i},3)
        X = [X;MDP.xn{i}(:,:,j,j)];
    end
    XX = [XX X];
end

subplot(3,1,1)
imagesc(1 - XX'), colormap gray
title('Beliefs - second level')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('Time')
ylabel('States')

XXX = [];
P  = [];
for t = 1:length(MDP.mdp)
    XX = [];
    
    for i = 1:numel(MDP.mdp(t).xn)
        X = [];
        for j = 1:size(MDP.mdp(t).xn{i},3)
            X = [X;MDP.mdp(t).xn{i}(:,:,j,j)];
        end
        XX = [XX X];
    end
    XXX = [XXX; XX];
    P = [P MDP.mdp(t).R];
end

subplot(3,1,2)
imagesc(1 - P), colormap gray
title('Beliefs - policy')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('Time')
ylabel('Policies')

subplot(3,1,3)
imagesc(1 - XXX'), colormap gray
title('Beliefs - first level')
set(gca,'xtick',[])
set(gca,'ytick',[])
xlabel('Time')
ylabel('States')