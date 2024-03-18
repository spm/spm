function MDP = DEM_drone_vision
% Demo of domain factors in the setting of active vision
%__________________________________________________________________________
%
% This routine illustrates weight sharing in the context of partially
% observed Markov decision processes. Specifically, it considers the
% problem of building posterior beliefs (c.f., a cognitive map) based upon
% spas sampling of an (unchanging) environment using the same likelihood
% mappings for each modality that conditioning the subset of hidden factors
% (i.e., the states of some location) that constitute the domain (i.e.,
% parents)  of the likelihood mapping. In other words, contextualising the
% allocentric representation of the environment to generate egocentric
% observations. Effectively, this means that one can encode the requisite
% projective geometry, not in the likelihood mappings, but in the hidden
% states to which they refer.
%
% The particular problem here models the environment as a collection of
% conditionally independent factors, where each factor corresponds to a
% location (a point on a 3D grid representation of the world). The levels
% of each factor encode the state or (semantic) class of each location. The
% conditional independence of the latent factors means that a subset of
% factors are the parents of the likelihood mapping at any one time. And
% these parents depend upon where a drone is currently looking.
% 
% A drone can fly around the environment with several lines of sight
% (corresponding to the number of modalities). Depending upon the location
% and orientation of the drone, a small number of locations are specified
% to generate what would be seen from that perspective. This means there
% are hidden factors pertaining to the location and orientation of the
% drone and a set of hidden factors corresponding to the state of the
% environment each location. The former are controllable, the latter are
% not. In this setup, the configuration of the environment is encoded by
% beliefs and exploration or navigation rests upon the expected information
% gain about states (as opposed to parameters). One can also include an
% expected cost of flying too close to locations with certain attributes.
%
% Note that this requires two sets of variational updates — to update
% posteriors over the state of the drone (based upon the likelihood of
% observations from different perspectives), and to update posteriors over
% the state of the environment (based upon the likelihood of observations,
% given that perspective). In principle, this allows the drone to infer
% both the environment and its relationship to that environment; although,
% in practice, we will illustrate the behaviour of a drone given precise
% beliefs about its initial state, or precise beliefs about the world
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $


%% set up and preliminaries
%==========================================================================
% First, specifying the size of the environment, the number and deployment
% of lines of sight and the number of directions the drone can turn among.
%==========================================================================
rng(1)
T   = 32;                                % number of moves
N   = 0;                                 % depth of planning
Nx  = 16;                                % size of environment
Ny  = 16;                                % size of environment
Nz  = 4;                                 % size of environment
Nd  = 4;                                 % depth of line of sight
Nc  = 4;                                 % number of classes
Nr  = 5;                                 % number of lines of sight
Na  = 8;                                 % number of drone angles

% Create a random environment where the class or state of each location
% increases with height
%--------------------------------------------------------------------------
S   = spm_conv(rand(Nx,Ny),2,2);
S   = S/max(S(:));
S   = round(S*(Nc + 2)) - 2;
S   = max(S,1);

for z = 1:Nz
    W(:,:,z) = (S == z)*z;
end

% boundary conditions (e.g., walls and floors; the last class)
%--------------------------------------------------------------------------
W    = max(W,1);
W(1,:,:)   = Nc;
W(end,:,:) = Nc;
W(:,1,:)   = Nc;
W(:,end,:) = Nc;
W(:,:,1)   = Nc;

% Now specify the hidden states at each location (a MAP)
%--------------------------------------------------------------------------
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            D{i,j,k} = sparse(W(i,j,k),1,1,Nc,1);
        end
    end
end
MAP  = D(:);
D    = MAP;

spm_figure('GetWin','Active inference'); clf
subplot(2,2,2), spm_show_x(MAP,Nx,Ny,Nz)
title('True environment'), drawnow

% assume a stationary world, in terms of transitions and initial states
%--------------------------------------------------------------------------
for f = 1:numel(MAP)
    B{f}  = spm_dir_norm(eye(Nc,Nc));
    D{f}  = spm_softmax(MAP{f}/8);
end

% Likelihood mapping: given states along the line of sight (i.e., rays)
%==========================================================================
% Given the state of locations at various depths along a line of
% sight, specifying the outcome in terms of depth and class; namely, the
% state of the first non-empty location, where the first state corresponds
% to empty space.
%--------------------------------------------------------------------------
Ns    = kron(ones(1,Nd),Nc);
depth = false([Nd + 1,Ns]);
state = false([Nc,    Ns]);
for s = 1:prod(Ns)

    % combination of classes
    %----------------------------------------------------------------------
    c   = spm_index(Ns,s);                 % states at  increasing depths
    d   = find(c > 1,1,'first');           % first nonempty state
    ind = num2cell(c);
    if numel(d)
        depth(d,   ind{:}) = true;         % depth of first occluder
        state(c(d),ind{:}) = true;         % state of first occluder
    else
        depth(end, ind{:}) = true;         % nothing in range
        state(1,   ind{:}) = true;         % nothing in range (empty)
    end
end

% add uncertainty about likelihood
%--------------------------------------------------------------------------
depth = depth + 1/8;
depth = spm_dir_norm(depth);
state = state + 1/8;
state = spm_dir_norm(state);

% lines of sight : the same for every modality or ray
%--------------------------------------------------------------------------
Ng    = Nr*Nr;
for g = 1:Ng
    A{g}      = depth;                    % first Ng modalities
    A{g + Ng} = state;                    % final Ng modalities
end

% Now specify the hidden states of the drone (and their dynamics)
%--------------------------------------------------------------------------
u     = [-1,0,1];
Nu    = numel(u);
for i = 1:Nu
    BX(:,:,i) = full(spm_speye(Nx,Nx,u(i),2));
    BY(:,:,i) = full(spm_speye(Ny,Ny,u(i),2));
    BZ(:,:,i) = full(spm_speye(Nz,Nz,u(i),2));
    BR(:,:,i) = full(spm_speye(Na,Na,u(i),1));
end

B = [BR, B];
B = [BZ, B];
B = [BY, B];
B = [BX, B];

% add uncertainty about state transitions
%--------------------------------------------------------------------------
c     = 1/16;
for f = 1:4
    b           = B{f};
    B{f}(:,:,1) = b(:,:,1) + c*b(:,:,2) + c*b(:,:,1)^2;
    B{f}(:,:,2) = b(:,:,2) + c*b(:,:,1) + c*b(:,:,3);
    B{f}(:,:,3) = b(:,:,3) + c*b(:,:,2) + c*b(:,:,3)^2;
    B{f}        = spm_dir_norm(B{f});
end

% intial states
%--------------------------------------------------------------------------
c = round(Nx/2);
D = [full(sparse(1,1,1,Na,1)); D];           % drone angles
D = [full(sparse(3,1,1,Nz,1)); D];           % location (z)
D = [full(sparse(c,1,1,Ny,1)); D];           % location (y)
D = [full(sparse(c,1,1,Nx,1)); D];           % location (x)

% number of factors and number of states per factor
%--------------------------------------------------------------------------
Nf    = numel(D);
for f = 1:Nf
    Ns(f) = numel(D{f});
end


% for each state of the drone evaluate which location factors are the
% parents of the outcomes (i.e., the domain of the likelihood mapping)
%--------------------------------------------------------------------------
% Effectively, this is where one defines the projective geometry, in terms
% of the mapping from a three space to an outcome modality (i.e., sensor).
% In other words, here, we specify which hidden states are responsible for
% generating the outcome in any given modality.
%
% In this example, we assume the field-of-view has the same azimuthal and
% polar angle that spans the angular interval between each angular position
% the drone can occupy. For example, if the drone can point in six
% directions, the field of view subtends an angle of 60°. Furthermore, we
% assume this field of view points downwards such that the drone can what
% is in front of it.
%--------------------------------------------------------------------------
theta = linspace(-pi/6,pi/6,Nr);                     % ray angles
phi   = linspace(0,2*pi*(1 - 1/Na),Na);              % drone angles

% Likelihood domains (i.e., parents id.fg) for all drone and ray angles
%--------------------------------------------------------------------------
Nu    = 4;
Ng    = Nr*Nr;
id.ff = (1:Nu);
id.fg = cell(Ng + Ng,1);
Nsf   = Ns(1:Nu);

% location of hidden factors in rectangular coordinates
%--------------------------------------------------------------------------
for s = 1:prod([Nx,Ny,Nz])
    XYZ(s,:) = spm_index([Nx,Ny,Nz],s);
end

% for every combination of (Nu) domain factors (i.e., id.ff)
%--------------------------------------------------------------------------
for s = 1:prod(Nsf)

    % get centre and angles of drone
    %----------------------------------------------------------------------
    c      = spm_index(Nsf,s);
    origin = c(1:3);
    angle  = phi(c(4));

    % get sampled locations for this modality
    %----------------------------------------------------------------------
    for g = 1:Ng

        a     = spm_index([Nr,Nr],g);      % ray (line of sight) angle
        ph    = angle - theta(a(2));       % azimuthal angle
        th    = theta(a(1)) - theta(1);    % polar angle

        % for each depth along ray (i.e., line of sight)
        %------------------------------------------------------------------
        for d = 1:Nd

            % find nearest location in the latent state space
            %--------------------------------------------------------------
            xyz    = spm_spherical2rectang(origin,ph,th,d);
            [~,j]  = min(sum(minus(XYZ,xyz).^2,2));
            loc(d) = j + Nu;

        end

        % save parents of A{g} - latent locations
        %------------------------------------------------------------------
        ind = num2cell(c);
        id.fg{g     }{ind{:}} = loc;
        id.fg{g + Ng}{ind{:}} = loc;
    end
end


%% priors: (cost) C
%==========================================================================
% Make proximate (depth = 1) objects aversive
%--------------------------------------------------------------------------
for g = 1:numel(A)
    C{g} = spm_dir_norm(ones(size(A{g},1),1));
end
for g = 1:Ng
    C{g} = spm_softmax(sparse(1,1,-2,Nd + 1,1));
end

% Supplement with constraints in latent state space
%--------------------------------------------------------------------------
% Constraints in a high dimensional setting can be specified in a small
% number of latent factors by supplying a matrix of constraints, where each
% column corresponds to a distinct constraint and the number of rows
% corresponds to the number of hidden factors. This constraint matrix
% contains the index of the costly state for each factor. If an index is
% zero, the constraint is taken to be independent of the corresponding
% factor; otherwise, this conditional factor has to have a high posterior
% over the next state, before the constraint is implemented.
%--------------------------------------------------------------------------
Nf      = numel(B);
id.cid  = zeros(Nf,1);
c       = 0;

% in this illustration, the drone cannot occupy non-empty states
%--------------------------------------------------------------------------
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            f     = sub2ind([Nx,Ny,Nz],i,j,k) + Nu;
            for s = 2:Nc                      % for nonempty states
                c           = c + 1;          % next constraint
                id.cid(1,c) = i;              % contrained factor
                id.cid(2,c) = j;              % contrained factor
                id.cid(3,c) = k;              % contrained factor
                id.cid(f,c) = s;              % conditional factor
            end
        end
    end
end


% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify allowable actions (with
% an action for each controllable state
%--------------------------------------------------------------------------
u     = 1:Nu;
U     = zeros(1,Nf);
U(u)  = 1;

% true states
%--------------------------------------------------------------------------
s        = zeros(Nf,1);
s((Nu + 1):end) = W(:);

% MDP Structure, specifying T epochs of active vision
%==========================================================================
mdp.T = T;                        % numer of moves
mdp.U = U;                        % controllable actions
mdp.A = A;                        % likelihood probabilities
mdp.B = B;                        % transition probabilities
mdp.C = C;                        % prior constraints
mdp.D = D;                        % initial priors
mdp.N = N;                        % policy depth

mdp.id = id;                      % domains
mdp.s  = s;

% Solve - an example of surveillance
%==========================================================================
% In the first task, the drone has to explore its environs and form beliefs
% about the scene at hand. These beliefs can then be queried; for example,
% counting sheep (green locations).
%--------------------------------------------------------------------------
      tic;
MDP = spm_MDP_VB_XXX(mdp);
t   = toc;

fprintf('Compute time: %i ms/update\n',round(1000*t/MDP.T))

% illustrate belief updating
%--------------------------------------------------------------------------
spm_figure('GetWin','Belief updating'); clf
spm_MDP_VB_trial(MDP);


% illustrate scene construction and exploration
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference');
spm_behaviour(MDP,Nx,Ny,Nz)

%--------------------------------------------------------------------------
fprintf('How many green things did you find?\n\n')
%--------------------------------------------------------------------------
s     = 3;
n     = 0;
for f = (Nu + 1):Nf
    if MDP.X{f}(s,end) > .95
        n = n + 1;
    end
end
%--------------------------------------------------------------------------
fprintf('I am fairly confident I found at least %i\n\n',n)
%--------------------------------------------------------------------------

return


% Solve - an example of search and rescue
%==========================================================================
% In the second task, the agent has to find a green thing (i.e., class 3)
% and maintain surveillance over it. To simulate this, we specify a
% preference class 3 in the class or attribute modalities
%--------------------------------------------------------------------------
s     = 3;
for g = (1:Ng) + Ng
    c    = sparse(s,1,8,Nc,1);
    C{g} = spm_softmax(c);
end

MDP   = mdp;
MDP.C = C;
MDP.T = 16;
MDP   = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Search and rescue');
spm_behaviour(MDP,Nx,Ny,Nz)

%--------------------------------------------------------------------------
fprintf('I''ve found a big green thing\n\n')
%--------------------------------------------------------------------------


% Solve - an example of wayfinding with priors
%==========================================================================
% In this example, we Illustrate pathfinding using inductive inference
% under precise prior beliefs about latent states (i.e., cognitve map).
%--------------------------------------------------------------------------

% Equip agents with precise posterior beliefs about initial (unchanging)
% latent states (e.g., providing it with a map)
%--------------------------------------------------------------------------
MDP   = mdp;
MDP.T = Nx;
MDP.D(1 + Nu:end) = MAP;

% Specify intial (home) and final hidden (goal) states (hid)
%--------------------------------------------------------------------------
s1         = XYZ(find(W == 1,1,'first'),:);
s2         = XYZ(find(W == 1,1,'last'),:);

MDP.D{1}   = full(sparse(s1(1),1,1,Nx,1));   % inital state
MDP.D{2}   = full(sparse(s1(2),1,1,Ny,1));   % inital state
MDP.D{3}   = full(sparse(s1(3),1,1,Nz,1));   % inital state
MDP.id.hid = s2(:);                          % final  state
MDP        = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Wayfinding');
spm_behaviour(MDP,Nx,Ny,Nz)

%--------------------------------------------------------------------------
fprintf('Drop-off completed\n\n')
%--------------------------------------------------------------------------

% Solve - an example of self orientation
%==========================================================================
% In this example we switch precise priors over initial location of the
% drone for precise priors about the environment to simulate self-location
% in a GPS denied setting.
%--------------------------------------------------------------------------
spm_figure('GetWin','Orientation');

% imprecise priors over location and orientation
%--------------------------------------------------------------------------
MDP   = mdp;
for f = 1:Nu
    MDP.D{f} = ones(size(D{f}));
end

% supply precise priors over latent states
%--------------------------------------------------------------------------
MDP.D(1 + Nu:end) = MAP;

% active inference
%--------------------------------------------------------------------------
MDP.T = 8;
MDP   = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Orientation');
spm_MDP_VB_trial(MDP);

spm_figure('GetWin','Wayfinding');
spm_behaviour(MDP,Nx,Ny,Nz)


return

% Subroutines
%==========================================================================
function xyz   = spm_spherical2rectang(origin,ph,th,d)
% Returns the Cartesian coordinates from spherical coordinates (and origin)
% origin   - [x,y,z]
% ph       - azimuthal angle (radians) (in x-y plane)
% th       - polar (radians)
% d        - radial distance
%
%--------------------------------------------------------------------------
x   = origin(1) + d*cos(ph)*cos(th);
y   = origin(2) + d*sin(ph)*cos(th);
z   = origin(3) - d*sin(th);
xyz = [x,y,z];

return


function RGB = spm_colour(O)
% subfunction: returns an RGB rendering of a multinomial distribution
%--------------------------------------------------------------------------
MAP = [0 0 0;
    1 0 1;
    0 1 0;
    0 0 1;
    1 1 0;
    0 1 1;
    1 1 0;
    0 1 1;
    1 0 1];
MAP = MAP(1:numel(O),:)';
RGB = min(MAP*O,1);

return

function spm_show_x(x,Nx,Ny,Nz)
% display posteriors
% FORMAT spm_show_x(x,Nx,Ny,Nz)
% x  - likelihood tensors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% illustrate images
%--------------------------------------------------------------------------
u     = 1 - 1/16;
for t = 1:size(x,2)
    D     = reshape(x(:,t),Nx,Ny,Nz);
    for i = 1:Nx
        for j = 1:Ny
            for k = 1:Nz

                % definitive class
                %----------------------------------------------------------
                [m,d] = max(D{i,j,k});
                if m > u
                    W(i,j,k) = d;
                else
                    W(i,j,k) = NaN;
                end
            end
        end
    end

    % create coloured image
    %----------------------------------------------------------------------
    I     = zeros(Nx,Ny,3);
    for i = 1:Nx
        for j = 1:Ny

            % find first definitive class
            %--------------------------------------------------------------
            s = squeeze(W(i,j,:));
            k = find(isfinite(s) & s > 1,1,'last');
            if numel(k)
                I(i,j,:) = spm_colour(D{i,j,k});
            end
        end
    end

    % display
    %----------------------------------------------------------------------
    imagesc(I), axis image

end

return

function spm_behaviour(MDP,Nx,Ny,Nz)
% display posteriors
% FORMAT spm_show_x(x,Nx,Ny)
% a{g} - likelihood tensors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% drone's location and observations
%--------------------------------------------------------------------------
[Nf,Ns,Nu,Ng,No] = spm_MDP_size(MDP);

Nr    = Ng/2;                     % number of rays
Nd    = sqrt(Nr);                 % number of rays (per dimension)
N     = prod([Nx,Ny,Nz]);         % number of locations
Nu    = Nf - N;                   % number of domain factors
X     = MDP.s(1:Nu,:);            % states
O     = MDP.o;                    % outcomes
L     = (1:N) + Nu;               % indices of locations

% locations
%--------------------------------------------------------------------------
for s = 1:prod([Nx,Ny,Nz])
    XYZ(s,:) = spm_index([Nx,Ny,Nz],s);
end

% illustrate behaviour
%--------------------------------------------------------------------------
for t = 1:MDP.T

    % true scene
    %----------------------------------------------------------------------
    Q     = {};
    for f = 1:Nf
        Q{f,1} = sparse(MDP.s(f,t),1,1,Ns(f),1);
    end
    subplot(2,2,2), hold off
    spm_show_x(Q(L),Nx,Ny,Nz), hold on
    title('True scene')

    % inferred scene
    %----------------------------------------------------------------------
    Q    = MDP.X;
    for f = 1:Nf
        Q{f} = Q{f}(:,t);
    end
    subplot(2,2,1), hold off
    spm_show_x(Q(L)',Nx,Ny,Nz), hold on
    title(sprintf('Inferred scene (t = %i)',t))

    % where the drone thinks it is
    %----------------------------------------------------------------------
    for f = 1:Nu
        [~,j] = max(Q{f});
        Xt(f) = j;
    end
    plot(Xt(2), Xt(1), '*y','MarkerSize',32,'LineWidth',4)

    % true location and orientation of drone
    %----------------------------------------------------------------------
    plot(X(2,t),X(1,t),'*w','MarkerSize',32,'LineWidth',4)

    % where the drone has been
    %----------------------------------------------------------------------
    plot(X(2,1:t),X(1,1:t),'.w','MarkerSize',16)

    % where the drone looks
    %----------------------------------------------------------------------
    j     = unique(spm_vec(MDP.j{:,t}));
    for i = 1:numel(j)
        x = spm_index([Nx,Ny,Nx],j(i) - Nu);
        plot(x(2),x(1),'ow','MarkerSize',2*x(3))
    end
   
    % 3D view
    %======================================================================
    subplot(4,2,6), hold off
    for s = 2:Ns(end)
        i = find(MDP.s((1:N) + Nu,t) == s);
        c = spm_colour(sparse(s,1,1,Ns(end),1));
        plot3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'o','MarkerSize',4,'Color',c)
        hold on
    end

    % add sampled locations
    %----------------------------------------------------------------------
    plot3(X(1,t),X(2,t),X(3,t),'*k','MarkerSize',16,'LineWidth',2), hold on
    for f = 1:numel(j)
        x = spm_index([Nx,Ny,Nx],j(f) - Nu);
        plot3(x(1),x(2),x(3),'ok','MarkerSize',2)
    end
    view(60,30), axis image
    title('Field of view')

    subplot(4,2,7), hold off
    plot(X(2,t),  X(3,t),  '*y','MarkerSize',16,'LineWidth',2), hold on
    plot(X(2,1:t),X(3,1:t),'.k','MarkerSize',16)
    axis image, axis([1 Ny 1 Nz]) 
    title('Height')

    % what the drone sees (depth)
    %----------------------------------------------------------------------
    g     = 1:Nr;
    depth = MDP.o(g,t);
    depth = reshape(depth,Nd,Nd);
    subplot(4,4,9), hold off
    imagesc(1 - depth), axis image
    title('Depth')

    % what the drone sees (state or semantic class)
    %----------------------------------------------------------------------
    g     = (1:Nr) + Nr;
    state = MDP.o(g,t);
    state = reshape(state,Nd,Nd);
    subplot(4,4,9), hold off

    % find class
    %--------------------------------------------------------------
    I     = zeros(Nd,Nd,3);
    for i = 1:Nd
        for j = 1:Nd
            D        = sparse(state(i,j),1,1,No(end),1);
            I(i,j,:) = spm_colour(D);
        end
    end
    subplot(4,4,10), hold off
    imagesc(I), axis image
    title('Class')

    % save movie
    %----------------------------------------------------------------------
    drawnow
    MOV(t) = getframe(gcf);

end

% Place movie in graphic subject
%--------------------------------------------------------------------------
set(gcf,'Userdata',{MOV,8})
set(gcf,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

return




