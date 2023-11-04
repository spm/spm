function MDP = DEM_drone
% Demo of structure learning (i.e.,disentanglement) of MNIST digits
%__________________________________________________________________________
%
% This routine illustrates weight sharing in the context of partially
% observed Markov decision processes. Specifically, it considers the
% problem of building posterior beliefs (c.f., a cognitive map) based upon
% spas sampling of an (unchanging) environment using the same likelihood
% mappings for each modality that conditioning the subset of hidden factors
% (i.e., the states of some location) that constitute the domain of the
% likelihood mapping. In other words, contextualising the allocentric
% representation of the environment to generate egocentric observations.
% Effectively, this means that one can encode the requisite projective
% geometry, not in the likelihood mappings, but in the hidden states to
% which they refer.
% 
% The particular problem illustrated here treats the environment as some
% surface; encoded by the height of a service at a number of grid points or
% locations. Each location corresponds to a hidden factor and the height is
% encoded by the state of the factor. A drone can fly around the
% environment with several lines of sight (corresponding to the number of
% modalities). Depending upon the location and orientation of the drone, a
% small number of locations are specified to generate what would be seen
% from that perspective. This means there are hidden factors pertaining to
% the location and orientation of the drone and a set of hidden factors
% corresponding to the state of the environment each location. The former
% are controllable, the latter are not. In this setup, the configuration of
% the environment is encoded by beliefs and exploration or navigation rests
% upon the expected information gain about states (as opposed to
% parameters). One can also include an expected cost of flying near high
% domains in the environment.
% 
% Note that this requires two sets of variational updates — to update
% posteriors over the state of the drone (based upon the likelihood of
% observations from different perspectives), and to update posteriors over
% the state of the environment (based upon the likelihood of observations,
% given that perspective). In principle, this allows the drone to infer
% both the environment and its relationship to that environment; although,
% in practice, we will illustrate the behaviour of a drone given precise
% beliefs about its initial state.
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $


%% set up and preliminaries
%==========================================================================
% First, specify the environment as a surface that can be at ground level,
% low or high (where high locations in the near field of view (i.e., line
% of sight) are, a priori, surprising or costly.
%==========================================================================
rng(1)
T   = 128;                               % number of moves
Nx  = 32;                                % size of environment
Ny  = 32;                                % size of environment
Nz  = 3;                                 % number of heights 
Nd  = 3;                                 % depth of line of sight

S   = spm_conv(rand(Nx,Ny),2,2);
S   = round(S/max(S(:))*Nz);
S   = max(S,1);

% boundary conditions (aversive: high/red)
%--------------------------------------------------------------------------
S(1,:)   = Nz;
S(end,:) = Nz;
S(:,end) = Nz;
S(:,1)   = Nz;

spm_figure('GetWin','Active inference'); clf
subplot(2,2,4), imagesc(S), title('True background'), axis image
subplot(2,2,2), bar3(S), title('True background'), axis image

% Now specify the hidden states at each location
%--------------------------------------------------------------------------
for i = 1:Nx
    for j = 1:Ny
        D{i,j} = sparse(S(i,j),1,1,Nz,1);
    end
end
D     = D(:);
subplot(2,2,3), spm_show_x(D,Nx,Ny), title('True background')

% stationary backgound
%--------------------------------------------------------------------------
for f = 1:numel(D)
    B{f} = eye(Nd,Nd);
    D{f} = spm_softmax(D{f}/8);
end


% Now specify the likelihood mapping from any location
%--------------------------------------------------------------------------
% Effectively, this is a specification of the multidimensional observations
% obtained in any given outcome modality, given the hidden states. In this
% example, each modality is a line of sight, and the outcomes are 10
% arrangements; for example, a tall (red) thing is in the distance (last
% sample point in the line of sight), while a short (green) thing is
% directly in front (i.e., the drone is hovering over a green thing). This
% would be outcome 5 below.
%--------------------------------------------------------------------------
a      = zeros([4,kron(ones(1,Nd),Nz)]);
for s1 = 1:size(a,2)
    for s2 = 1:size(a,3)
        for s3 = 1:size(a,4)

            % tall location
            %--------------------------------------------------------------
            if s1 == Nd
                o = 1;             % tall-front
            elseif s2 == Nd
                if s1 == 1
                    o = 2;         % tall-near
                elseif s1 > 1
                    o = 3;         % tall-near short-front
                end
            elseif s3 == Nd
                if s1 == 1
                    if s2 == 1
                        o = 4;     % tall-back
                    else
                        o = 5;     % tall-back short-near
                    end
                elseif s1 > 1
                    o = 6;         % tall-back short-front
                end

                % no tall location
                %----------------------------------------------------------
            else
                if s1 == 1
                    if s2 == 1
                        if s3 == 1
                            o = 7; % nothing
                        else
                            o = 8; % short-back
                        end
                    else
                        o = 9;     % short-near
                    end
                else
                    o = 10;        % short-front
                end
            end

            % liklihood tensor
            %--------------------------------------------------------------
            a(o,s1,s2,s3) = 1;

        end
    end
end

% three lines of sight
%--------------------------------------------------------------------------
Ng    = 3;
for g = 1:Ng
    A{g} = a;
end

% Now specify the hidden states of the drone (and their dynamics)
%--------------------------------------------------------------------------
c     = round(Nx/2);
Nr    = 4;
u     = [-1,0,1];
Nu    = numel(u);
for i = 1:Nu
   BX(:,:,i) = full(spm_speye(Nx,Nx,u(i),2));
   BY(:,:,i) = full(spm_speye(Ny,Ny,u(i),2));
   BR(:,:,i) = full(spm_speye(Nr,Nr,u(i),1));
end

B = [BR, B];
B = [BY, B];
B = [BX, B];
D = [full(sparse(1,1,1,Nr,1)); D];
D = [full(sparse(c,1,1,Ny,1)); D];
D = [full(sparse(c,1,1,Nx,1)); D];


% for each state of the drone evaluate which location factors are the
% domain of the likelihood mapping
%--------------------------------------------------------------------------
% Effectively, this is where one defines the perspective geometry, in terms
% of the mapping from a three space to an outcome modality (i.e., sensor).
% In this example, the drone has three lines of sight starting beneath it
% and progressing on the ground plane in three directions (directly ahead,
% 45° to the drone’s left and to the drone’s right). In other words, here,
% we specify which hidden states are responsible for generating the outcome
% in any given modality.
%--------------------------------------------------------------------------
Nf    = numel(B);
Ng    = 3;
fg    = cell(Ng,1);
for x = 1:Nx
    for y = 1:Ny

        xp = min( (1:Nd) + x - 1,Nx);
        yp = min( (1:Nd) + y - 1,Ny);
        xn = max( 1 - (1:Nd) + x,1);
        yn = max( 1 - (1:Nd) + y,1);
        x0 = zeros(1,Nd) + x;
        y0 = zeros(1,Nd) + y;

        % North pointing
        %------------------------------------------------------------------
        fg{1}{x,y,1} = sub2ind([Nx,Ny],xp,y0) + Nu;
        fg{2}{x,y,1} = sub2ind([Nx,Ny],xp,yp) + Nu;
        fg{3}{x,y,1} = sub2ind([Nx,Ny],x0,yp) + Nu;

        % East pointing
        %------------------------------------------------------------------
        fg{1}{x,y,2} = sub2ind([Nx,Ny],x0,yp) + Nu;
        fg{2}{x,y,2} = sub2ind([Nx,Ny],xn,yp) + Nu;
        fg{3}{x,y,2} = sub2ind([Nx,Ny],xn,y0) + Nu;

        % South pointing
        %------------------------------------------------------------------
        fg{1}{x,y,3} = sub2ind([Nx,Ny],xn,y0) + Nu;
        fg{2}{x,y,3} = sub2ind([Nx,Ny],xn,yn) + Nu;
        fg{3}{x,y,3} = sub2ind([Nx,Ny],x0,yn) + Nu;

        % West pointing
        %------------------------------------------------------------------
        fg{1}{x,y,4} = sub2ind([Nx,Ny],x0,yn) + Nu;
        fg{2}{x,y,4} = sub2ind([Nx,Ny],xp,yn) + Nu;
        fg{3}{x,y,4} = sub2ind([Nx,Ny],xp,y0) + Nu;

    end
end

% domains
%--------------------------------------------------------------------------
u     = (1:Nu);
id.fg = fg;
id.ff = u;


%% priors: (cost) C
%==========================================================================
% Finally, we have to specify the prior constraints in terms of
% probabilities over outcomes: first aversive (red) outcome is tall-front
%--------------------------------------------------------------------------
c     = zeros(size(a,1),1);
c(1)  = -32;                        
for g = 1:Ng
    C{g} = spm_softmax(c);
end

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify allowable actions (with
% an action for each controllable state
%--------------------------------------------------------------------------
U     = zeros(1,Nf);
U(u)  = 1;

% true states
%--------------------------------------------------------------------------
s      = zeros(Nf,T);
i      = (1:(Nx*Ny)) + Nu;
s(i,:) = S(:)*ones(1,T);


% MDP Structure, specifying 8 epochs (i.e., 4 seconds of active vision)
%==========================================================================
mdp.T = T;                        % numer of moves
mdp.U = U;                        % controllable actions
mdp.A = A;                        % likelihood probabilities
mdp.B = B;                        % transition probabilities
mdp.C = C;                        % prior constraints
mdp.D = D;                        % initial priors
mdp.N = 0;                        % policy depth

mdp.id = id;                      % domains
mdp.s  = s;

% Solve - an example with multiple epochs
%==========================================================================
MDP    = spm_MDP_VB_XXX(mdp);

% illustrate belief updating of second agent
%--------------------------------------------------------------------------
spm_figure('GetWin','Beleif updating'); clf
spm_MDP_VB_trial(MDP);


% illustrate scene construction and perceptual synthesis
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference');
subplot(2,2,1)
spm_behaviour(MDP,Nx,Ny)


return

% Subroutines
%==========================================================================

function spm_show_x(x,Nx,Ny)
% display posteriors
% FORMAT spm_show_x(x,Nx,Ny)
% a{g} - likelihood tensors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% illustrate images
%--------------------------------------------------------------------------
for t = 1:size(x,2)
    D     = reshape(x(:,t),Nx,Ny);
    for i = 1:Nx
        for j = 1:Ny
            I(i,j,:) = flip(full(D{i,j}));
        end
    end
    imagesc(I), axis image
    drawnow
end

return

function spm_behaviour(MDP,Nx,Ny)
% display posteriors
% FORMAT spm_show_x(x,Nx,Ny)
% a{g} - likelihood tensors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% drone's location and observations
%--------------------------------------------------------------------------
Nu    = 3;
X     = MDP.s(1:Nu,:);
O     = MDP.o;
L     = (1:(Nx*Ny)) + Nu;

% illustrate behaviour
%--------------------------------------------------------------------------
for t = 1:MDP.T

    % Inferred background
    %----------------------------------------------------------------------
    Q    = MDP.X;
    for q = 1:numel(Q)
        Q{q} = Q{q}(:,t);
    end
    cla, spm_show_x(Q(L)',Nx,Ny), hold on

    % Location and orientation of drone
    %----------------------------------------------------------------------
    d     = 1/2;
    plot(X(2,t),X(1,t),'*k','MarkerSize',32)
    if X(3,t) == 1, plot(X(2,t) + 0,X(1,t) + d,'ok','MarkerSize',16), end
    if X(3,t) == 2, plot(X(2,t) + d,X(1,t) + 0,'ok','MarkerSize',16), end
    if X(3,t) == 3, plot(X(2,t) - 0,X(1,t) - d,'ok','MarkerSize',16), end
    if X(3,t) == 4, plot(X(2,t) - d,X(1,t) - 0,'ok','MarkerSize',16), end

    % where the drone has been
    %----------------------------------------------------------------------
    plot(X(2,1:t),X(1,1:t),'.w','MarkerSize',16)

    % where the drone looks
    %----------------------------------------------------------------------
    j     = unique(spm_vec(MDP.j{:,t}));
    [x,y] = ind2sub([Nx,Ny],j - Nu);
    plot(y,x,'ow','MarkerSize',8)

    title('Inferred background')
    drawnow
    drone(t) = getframe(gca);

end

set(gca,'Userdata',{drone,8})
set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

return




