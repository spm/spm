function MDP = DEM_drone_inspection
% Demo of domain factors in the setting of active vision
%__________________________________________________________________________
%
% This demonstration routine illustrates active inference in the service of
% inspecting objects in the scene; specifically, furniture or plant in a
% room. The routine demonstrates three aspects of the ensuing behaviour.
% First, the agent (e.g., cleaning robot) explores the room to form a
% cognitive map of where things are. Objects (e.g., tables) that are
% cleanable are recognised in terms of belonging to a cleanable class that
% may be clean (dark green) or unclean (light green). Having established
% precise posterior beliefs about the scene, the agent is equipped with
% prior beliefs that each cleanable location (i.e., surface) can either be
% clean or unclean and is given prior preferences for seeing (and
% implicitly cleaning or spraying) unclean surfaces. The agent then
% resolves uncertainty about whether surfaces are clean or unclean and
% observes them to be clean when and only when, they have been observed
% (and implicitly cleaned). This is because the state of the world (i.e.,
% generative process) is always clean; however, the agent is not convinced
% about this until it has seen clean surfaces for itself. Finally, the
% agent goes home after completing its cleaning task using inductive
% inference to specify its sleeping place (avoiding obstacles in the usual
% way).
%
% This example assumes local telemetry is available; in other words, the
% position relative to some arbitrary (allocentric) frame of reference is
% known — and the agent is equipped with both depth and RGB sensors, where
% the RGB data are (semantically) segmented into eight classes using some
% amortisation (computer vision) scheme. The agent has four degrees of
% freedom. It can move on the floor in two dimensions. It%s sensors are
% mounted on a vertical pole and can be moved up and down and rotated.
%
%__________________________________________________________________________
% Copyright (C) 2019 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_surveillance.m 8313 2022-09-30 18:33:43Z karl $


%% set up and preliminaries
%==========================================================================
% First, specify the size of the environment, the number and deployment of
% lines of sight and the number of directions sensors can point.
%==========================================================================
rng(1)

Nx  = 32;                                % size of environment
Ny  = 28;                                % size of environment
Nz  = 8;                                 % size of environment
Nd  = 16;                                % depth of rays
Nc  = 8;                                 % number of classes
Nr  = 16;                                % number of rays
Na  = 12;                                % number of camera angles

% global parameters
%--------------------------------------------------------------------------
global PG %#ok<GVMIS> 
PG.Nx  = Nx;
PG.Ny  = Ny;
PG.Nz  = Nz;
PG.Nd  = Nd;
PG.Nr  = Nr;
PG.Nc  = Nc;
PG.Na  = Na;

PG.phi = linspace(0,2*pi*(1 - 1/Na),Na); % camera angles
PG.FOV = linspace(-pi/6,pi/6,Nr);        % ray angles (60 degree FOV)

% Get environment
%--------------------------------------------------------------------------
[D,W]  = spm_drone_unity(Nx,Ny,Nz,Nc);
MAP    = D(:);
D      = MAP;

% 2D view: true scene
%==========================================================================
spm_figure('GetWin','Active inference'); clf
subplot(2,2,2), spm_show_x(MAP,Nx,Ny,Nz)
title('True environment')

% 3D view: true scene
%==========================================================================

% location of hidden factors in rectangular coordinates
%--------------------------------------------------------------------------
XYZ   = spm_combinations([Nx,Ny,Nz]);
subplot(4,2,6), hold off
for s = 2:Nc
    i = find(W(:) == s);
    c = spm_colour(sparse(s,1,1,Nc,1));
    plot3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'.','MarkerSize',32,'Color',c)
    hold on
end
axis image
drawnow

% Likelihood mapping: given states along the line of sight (i.e., rays)
%==========================================================================
% Given the state of locations at various depths along a line of sight,
% specify the outcome in terms of depth and class; namely, the state of the
% first non-empty location, where the first state corresponds to empty
% space. The requisite likelihood tensors are specified with function
% handles (see subroutines)
%--------------------------------------------------------------------------

% lines of sight: the same for every modality or ray: process
%--------------------------------------------------------------------------
Ng    = Nr*Nr;
for g = 1:Ng
    GA{g}       = @spm_depth;              % 1st Ng modalities
    GA{g + Ng}  = @spm_state;              % 2nd Ng modalities
end

% add telemetry (controllable modalities)
%--------------------------------------------------------------------------
GA{Ng + Ng + 1} = eye(Nx,Nx);              % x - location (x)
GA{Ng + Ng + 2} = eye(Ny,Ny);              % y - location (y)
GA{Ng + Ng + 3} = eye(Nz,Nz);              % z - location (z)
GA{Ng + Ng + 4} = eye(Na,Na);              % a - camera angles

% controlled outcomes (telemetry)
%--------------------------------------------------------------------------
ID.control = (1:4) + Ng + Ng;

% for every combination of (Nu) domain factors (i.e., id.ff) specify
% projective geometry
%--------------------------------------------------------------------------
Nu    = 4;
ID.ff = 1:Nu;
ID.fg = @spm_fg;

% repeat for generative model
%==========================================================================

% depth and vision modalities: model
%--------------------------------------------------------------------------
Ng    = Nr*Nr;
A     = cell(Ng + Nu,1);
for g = 1:Ng
    A{g}      = @A_depth;                 % 1st Ng modalities
    A{g + Ng} = @A_state;                 % 2nd Ng modalities
end

% add telemetry
%--------------------------------------------------------------------------
A{Ng + Ng + 1} = eye(Nx,Nx);              % x - location (x)
A{Ng + Ng + 2} = eye(Ny,Ny);              % y - location (y)
A{Ng + Ng + 3} = eye(Nz,Nz);              % z - location (z)
A{Ng + Ng + 4} = eye(Na,Na);              % a - camera angles

% indices of central and peripheral (visual) modalities for planning
%--------------------------------------------------------------------------
xy    = spm_combinations([Nr,Nr]);
ge    = find(sum(minus(xy(:,2),Nr/2).^2,2) < 1);

% modalities subtending EFE
%--------------------------------------------------------------------------
id.ge = unique([ge; (Ng + ge)]');

% uncomment to remove depth modality
%--------------------------------------------------------------------------
% id.g = {Ng:(Ng + Ng + 4)};

% convert likelihood arrays to logical arrays where possible
%--------------------------------------------------------------------------
for g = 1:numel(GA)
    if isnumeric(GA{g})
        GA{g} = logical(GA{g});
    end
end
for g = 1:numel(A)
    if isnumeric(GA{g})
        A{g} = logical(A{g});
    end
end

% latent states:
%==========================================================================
% x - location (x)
% y - location (y)
% z - location (z)
% a - camera angle
% s - state at location 1
% ...
% s - state at location Nx*Ny*Nz = Nf
%--------------------------------------------------------------------------
Nu    = 4;                                 % # controllable factors
id.ff = 1:Nu;                              % domain factors (x,y,x,a)
id.fg = @spm_fg;                           % parents of A{g}

% Specify transition priors (i.e., dynamics)
%--------------------------------------------------------------------------
d     = 4;
u     = [-1,0,1];
for i = 1:numel(u)

   % perimeter contraints on x
    %----------------------------------------------------------------------
    BX(:,:,i) = full(spm_speye(Nx,Nx,u(i),1));
    for j = 1:Nx
        if j < d
            BX(:,j,i) = full(sparse(d,1,1,Nx,1));
        end
        if j > (Nx - d + 1)
            BX(:,j,i) = full(sparse((Nx - d + 1),1,1,Nx,1));
        end
    end

    % perimeter contraints on y
    %----------------------------------------------------------------------
    BY(:,:,i) = full(spm_speye(Ny,Ny,u(i),2));
    for j = 1:Ny
        if j < d
            BY(:,j,i) = full(sparse(d,1,1,Ny,1));
        end
        if j > (Ny - d + 1)
            BY(:,j,i) = full(sparse((Ny - d + 1),1,1,Ny,1));
        end
    end

    % height contraints on z
    %----------------------------------------------------------------------
    BZ(:,:,i) = full(spm_speye(Nz,Nz,u(i),2));
    for j = 1:Nz
        if j < 3
            BZ(:,j,i) = full(sparse(3,1,1,Nz,1));
        end
        if j > (Nz - 1)
            BZ(:,j,i) = full(sparse((Nz - 1),1,1,Nz,1));
        end
    end

    % rotational moves
    %----------------------------------------------------------------------
    BA(:,:,i) = full(spm_speye(Na,Na,u(i),1));

end

% assume a stationary world, in terms of transitions and initial states
%--------------------------------------------------------------------------
GB    = {BX BY BZ BA};                           % controllable factors
for f = 1:(Nx*Ny*Nz)                             % latent place factors
    GB{end + 1}  = logical(eye(Nc,Nc));
end

% add uncertainty about state transitions (noise)
%--------------------------------------------------------------------------
n     = [1/16, 1/16, 1/64, 0]*exp(-16);          % state noise
for f = 1:numel(n)
    b = GB{f};
    c = n(f);
    GB{f}(:,:,1) = b(:,:,1) + c*b(:,:,2) + c*b(:,:,1)^2;
    GB{f}(:,:,2) = b(:,:,2) + c*b(:,:,1) + c*b(:,:,3);
    GB{f}(:,:,3) = b(:,:,3) + c*b(:,:,2) + c*b(:,:,3)^2;
    GB{f}        = spm_dir_norm(GB{f});
end

% priors over initial states
%==========================================================================
B     = GB;                                      % state transitions
for f = 1:numel(B)
    D{f} = ones(size(B{f},1),1);
end

V    = fix([Nx - 4,Ny - 4,Nz - 2,Na/2]);         % initial location
D{1} = full(sparse(V(1),1,1,Nx,1));              % location (x)
D{2} = full(sparse(V(2),1,1,Ny,1));              % location (y)
D{3} = full(sparse(V(3),1,1,Nz,1));              % location (z)
D{4} = full(sparse(V(4),1,1,Na,1));              % camera angles


% priors: (cost) C
%==========================================================================

% no modality contraints at this stage
%--------------------------------------------------------------------------
for g = 1:numel(A)
    C{g} = [];
end

% In addition, specify constraints in latent state-space (cif)
%--------------------------------------------------------------------------
id.cid = @spm_cid;

% This concludes the ABC of the model; namely, the likelihood mapping,
% prior transitions and preferences. Now, specify allowable actions (with
% an action for each controllable state
%--------------------------------------------------------------------------
u     = 1:Nu;
U     = zeros(1,numel(B));
U(u)  = 1;

% true states
%--------------------------------------------------------------------------
s     = [V(:); W(:)];

% MDP Structure, specifying T epochs of active vision
%==========================================================================
mdp.T  = 128;                     % number of moves
mdp.A  = A;                       % likelihood probabilities
mdp.B  = B;                       % transition probabilities
mdp.C  = C;                       % prior constraints
mdp.D  = D;                       % initial priors
mdp.N  = 0;                       % policy depth
mdp.U  = U;                       % controllable actions

mdp.GA = GA;                      % likelihood probabilities
mdp.GB = GB;                      % transition probabilities
mdp.GU = U;                       % controllable actions

mdp.ID = ID;                      % process domains
mdp.id = id;                      % model domains
mdp.s  = s;                       % inital states

% Solve - an example of surveillance
%==========================================================================
% In the first task, the agent has to explore its environs and form beliefs
% about the scene at hand.
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
spm_behaviour(MDP)


% Solve - an example of inspection
%==========================================================================
% In the second task, the agent has to search for unclean (class 3)
% services — and confirm that they are clean (class 4). To do
% this, we equip the agent with prior beliefs that cleanable locations
% could be clean or unclean and that it prefers to see unclean locations in
% its central field of view.
%--------------------------------------------------------------------------
c     = spm_softmax(sparse([3,4],1,[32,2],Nc,1));
for g = (1:Ng) + Ng
    C{g} = c;
end

% uninspected priors (UMAP)
%--------------------------------------------------------------------------
% (1) black (void)
% (2) blue  (floor)
% (3) green (unseen/unclean)
% (4) green (seen/clean)
% (5) brown (uncleanable)
% (6) pink  (wall)
% (7) brick (pillar)
% (8) red   (target)
%--------------------------------------------------------------------------
Nf    = prod([Nx,Ny,Nz]);
D     = MAP;
for f = 1:Nf
    if MAP{f}(4)
        D{f}(3) = 1;
        D{f}    = spm_dir_norm(D{f});
    end
end

% Equip agents with precise posterior beliefs about initial states
%--------------------------------------------------------------------------
f        = (1:Nf) + Nu;
MDP      = mdp;
MDP.T    = 32;
MDP.D(f) = D;
MDP.C    = C;
MDP      = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Inspection');
spm_behaviour(MDP)

% Solve - an example of wayfinding with priors
%==========================================================================
% In this example, we Illustrate pathfinding using inductive inference
% under precise prior beliefs about latent states (i.e., cognitve map).
%--------------------------------------------------------------------------

% Equip agents with precise posterior beliefs about initial states
%--------------------------------------------------------------------------
MDP      = mdp;
MDP.T    = 32;
MDP.D(f) = MAP;

% Specify intial (home) and final hidden (goal) states (hid)
%--------------------------------------------------------------------------
s0         = [16, 5, 4];                     % inital state
sT         = V;                              % final  state

MDP.s(1:4) = [s0(:); 1];
MDP.D{1}   = full(sparse(s0(1),1,1,Nx,1));   % inital state
MDP.D{2}   = full(sparse(s0(2),1,1,Ny,1));   % inital state
MDP.D{3}   = full(sparse(s0(3),1,1,Nz,1));   % inital state
MDP.D{4}   = full(sparse(    1,1,1,Na,1));   % inital state

MDP.id.hid = sT(:);                          % final  state
MDP        = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Go Home'); clf
spm_behaviour(MDP)
subplot(2,2,1), hold on, plot(sT(2),sT(1),'ow','MarkerSize',48)

return


% Subroutines
%==========================================================================

% Likelihood mapping: given states along the line of sight (i.e., rays)
%==========================================================================
% Given the state of locations at various depths along a line of
% sight, specify the outcome in terms of depth and class; namely, the
% state of the first non-empty location, where the first state corresponds
% to empty space.
%--------------------------------------------------------------------------

function O = spm_depth(s)
% likelihood of depth given some states (s)
% FORMAT O = spm_depth(s)
% s  - class or state along ray
%--------------------------------------------------------------------------
global PG

% depth outcome
%--------------------------------------------------------------------------
O   = zeros(PG.Nd + 1,1);
d   = find(s > 1,1,'first');           % first nonempty class
if numel(d)
    O(d)         = 1;                  % depth of first occluder
else
    O(PG.Nd + 1) = 1;                  % nothing in range
end
return


function O = spm_state(s)
% likelihood of class given some states (s)
% FORMAT O = spm_state(s)
% s  - class or state along ray
%--------------------------------------------------------------------------
global PG

% class (state)
%--------------------------------------------------------------------------
O   = zeros(PG.Nc,1);
d   = find(s > 1,1,'first');           % first nonempty class
if numel(d)
    O(s(d)) = 1;                       % class of first occluder
else
    O(1)    = 1;                       % nothing in range (null)
end
return


function [O,r] = A_depth(P,f,Q)
% likelihood of depth
% FORMAT [O,r] = A_depth(P,f,Q)
% P  - marginal over states along ray / outcomes
% f  - indices of nonzero elements
% Q  - marginal over states along ray
%__________________________________________________________________________
% This is a dual use function. If the first argument (P) is a cell, that
% this function returns the posterior predictive densities over outcomes,
% given that the posterior over states in P. Conversely, if P is a vector,
% then the marginal likelihoods over states are returned, given the
% outcomes in P. These marginals can be selected on the basis of posteriors
% over states in Q. The indices of the reduced marginals are returned in r.
%--------------------------------------------------------------------------
global PG

% depth outcome
%==========================================================================
if iscell(P)

    % outcome likelihood: O = spm_dot(A,P(n));
    %----------------------------------------------------------------------
    q = spm_cat(P);                        % marginals along ray
    q = q(1,:);                            % P(empty)
    p = 1 - q;                             % P(occupied)
    O = abs(p.*cumprod([1 q(1:end - 1)])); % P(first occupied)
    O = [O(:); 1 - sum(O)];                % P(outcome)

else

    % un-normalised likelihood of states: P = spm_dot(A,O);
    %======================================================================

    % return marginal likelihoods as cells
    %----------------------------------------------------------------------
    d  = find(P,1,'first');                % depth on this ray
    if d > PG.Nd
        o = false(PG.Nc,PG.Nd);
        o(1,:) = true;                     % all empty
    else
        o = false(PG.Nc,d);
        o(1,1:end - 1) = true;             % all empty
        o(2:end,end)   = true;             % until first occluder
    end

    % convert to cell array
    %----------------------------------------------------------------------
    d  = min(size(o,2),numel(f));
    O  = cell(1,d);
    for j = 1:d
        O{j} = o(f{j},j);
    end
    r  = 1:d;

end

return

function [O,r] = A_state(P,f,Q)
% likelihood of states
% FORMAT [O,r] = A_state(P,f,Q)
% P  - marginal over states along ray / outcomes
% f  - indices of nonzero elements
% Q  - marginal over states along ray
%--------------------------------------------------------------------------
% This is a dual use function. If the first argument (P) is a cell, that
% this function returns the posterior predictive densities over outcomes,
% given that the posterior over states in P. Conversely, if P is a vector,
% then the marginal likelihoods over states are returned, given the
% outcomes in P. These marginals can be selected on the basis of posteriors
% over states in Q. The indices of the reduced marginals are returned in r.
% depth outcome
%==========================================================================
global PG

if iscell(P)

    % outcome likelihood requested: O = spm_dot(A,P(n));
    %----------------------------------------------------------------------
    d            = A_depth(P);            % depth of first occluder
    q            = spm_cat(P);            % P(state)
    q(1,:)       = 0;                     % P(state | occupied)
    q(1,end + 1) = 1;                     % out of range

    q = spm_dir_norm(q);                  % normalise
    O = q*d;                              % P(outcome)

else

    % un-normalised likelihood of states: P = spm_dot(A,O);
    %----------------------------------------------------------------------
    % Only return a marginal likelihood for any state for which there is a
    % precise posterior belief it is at an occluding location: i.e., the
    % first occluder
    %----------------------------------------------------------------------
    O = {};
    r = [];
    d = A_depth(Q);                       % depth of first occluder
    d = find(d > (1 - 1/32),1,'first');
    if numel(d)
        if d <= PG.Nd && d <= numel(f)
            r = d;
            O = {P(f{r})};
        end
    end

end

return

function [D,hif] = spm_cid(Q)
% constraint function
% FORMAT [D,hif] = spm_cid(Q)
% D    - tensor of allowed latent states
% hif  - in subspace of hif factors
%--------------------------------------------------------------------------
global PG

% get empty locations
%--------------------------------------------------------------------------
Nu  = 4;
Nf  = PG.Nx*PG.Ny*PG.Nz;
Q   = spm_cat(Q((1:Nf) + Nu));

% allow if empty
%--------------------------------------------------------------------------
D   = Q(1,:) > (1 - 1/32);
D   = reshape(full(D),PG.Nx,PG.Ny,PG.Nz);
D   = spm_cross(all(D(:,:,(2:PG.Nz - 2)),3),ones(1,PG.Nz));

% (location) factors
%--------------------------------------------------------------------------
hif = 1:3;

return

function par = spm_fg(g,s)
% projective geometry for process
% FORMAT loc = spm_fg(g,s)
%
% s(1) - location (x)
% s(2) - location (y)
% s(3) - location (z)
% s(4) - camera angles
%
% par  - parents of g
%--------------------------------------------------------------------------
% Effectively, this is where one defines the projective geometry, in terms
% of the mapping from a three space to an outcome modality (i.e., sensor).
% In other words, here, we specify which hidden states are responsible for
% generating the outcome in any given modality.
%
% In this example, we assume the field-of-view has the same azimuthal and
% polar angle that spans the angular interval between each angular position
% the camera can occupy. For example, if the camera can point in six
% directions, the field of view subtends an angle of 60°. Furthermore, we
% assume this field of view points downwards such that the agent can what
% is in front of it.
%--------------------------------------------------------------------------
% for g = 1:Ng
%     GA{g}           = depth;               % 1st Ng modalities
%     GA{g + Ng}      = state;               % 2nd Ng modalities
% end
%
% telemetry
%--------------------------------------------------------------------------
% GA{Ng + Ng + 1} = eye(Nx,Nx);              % x - location (x)
% GA{Ng + Ng + 2} = eye(Ny,Ny);              % y - location (y)
% GA{Ng + Ng + 3} = eye(Nz,Nz);              % z - location (z)
% GA{Ng + Ng + 4} = eye(Na,Na);              % a - sensor angles
%--------------------------------------------------------------------------
global PG

Nu    = 4;                                   % domain factors
Ng    = PG.Nr*PG.Nr;                         % number of rays

% check modality
%--------------------------------------------------------------------------
if g > (Ng + Ng)

    % telemetry
    %----------------------------------------------------------------------
    par = g - (Ng + Ng);
    par = uint16(par);
    return

elseif g > Ng

    % class (same factors as depth)
    %----------------------------------------------------------------------
    g = g - Ng;

end

% factors on ray g: depth and class
%==========================================================================

% get centre and angles of sensor
%--------------------------------------------------------------------------
origin = [s(1:3)];
angle  = PG.phi(s(4));

% get sampled locations for this modality
%--------------------------------------------------------------------------
theta  = PG.FOV;                            % ray angles (60 degree FOV)
d      = spm_index([PG.Nr,PG.Nr],g);        % ray (line of sight) angle
ph     = angle - theta(d(2));               % azimuthal angle
th     = theta(d(1)) - theta(1);            % polar angle

% for each depth along ray (i.e., line of sight)
%--------------------------------------------------------------------------
par   = zeros(1,0);
for d = 1:PG.Nd

    % find nearest location in the latent state space
    %----------------------------------------------------------------------
    xyz = spm_spherical2rectang(origin,ph,th,d);
    x   = uint8(xyz(1));
    y   = uint8(xyz(2));
    z   = uint8(xyz(3));
    if x > 0 && x <= PG.Nx && y > 0 && y <= PG.Ny && z > 0 && z <= PG.Nz
        j      = spm_sub2ind([PG.Nx,PG.Ny,PG.Nz],x,y,z);
        par(d) = j + Nu;
    else
        break
    end
end

% save parents of A{g} - latent locations
%--------------------------------------------------------------------------
par   = uint16(par);

return

function xyz   = spm_spherical2rectang(origin,ph,th,d)
% Returns the Cartesian coordinates from spherical coordinates (and origin)
% origin   - [x,y,z]
% ph       - azimuthal angle (radians) (in x-y plane)
% th       - polar (radians)
% d        - radial distance
%--------------------------------------------------------------------------
x   = origin(1) + d*cos(ph)*cos(th);
y   = origin(2) + d*sin(ph)*cos(th);
z   = origin(3) - d*sin(th);
xyz = [x,y,z];

return


function RGB = spm_colour(O)
% subfunction: returns an RGB rendering of a multinomial distribution
%--------------------------------------------------------------------------
c   = 1/3;
MAP = [0 0 0;                     % (1) black
    169 196 219;                  % (2) blue
    027 101 073;                  % (3) green
    090 205 098;                  % (4) green
    204 204 000;                  % (5) brown
    255 229 204;                  % (6) pink
    204 102 000;                  % (7) brick
    255 000 000]/255;             % (8) red
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
u     = 2/numel(x{1});
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
    I     = ones(Nx,Ny,3)/2;
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

function [D,W] = spm_drone_unity(Nx,Ny,Nz,Nc)
% creates an indoor scene based on size of the environment
%==========================================================================

% Create a random environment with several objects of different classes
%--------------------------------------------------------------------------
% (1) black (void)
% (2) blue  (floor)
% (3) green (unseen/unclean)
% (4) green (seen/clean)
% (5) brown (uncleanable)
% (6) pink  (wall)
% (7) brick (pillar)
% (8) red   (target)

xyz      = spm_combinations([Nx,Ny,Nz]);
W        = ones(Nx,Ny,Nz);                % empty space
W(:,:,1) = 2;                             % floor / carpet

% place cleanable object in scene
%--------------------------------------------------------------------------
d    = 1 & (xyz(:,1) > 12) & (xyz(:,1) < (Nx - 4));
d    = d & (xyz(:,2) > 12) & (xyz(:,2) < (12 + 6));
d    = d & (xyz(:,3) > 4)  & (xyz(:,3) < 6);
W(d) = 4;

% place cleanable object in scene
%--------------------------------------------------------------------------
d    = 1 & (xyz(:,1) > 12) & (xyz(:,1) < (Nx - 4));
d    = d & (xyz(:,2) > 13) & (xyz(:,2) < (13 + 2));
d    = d & (xyz(:,3) > 0)  & (xyz(:,3) < 6);
W(d) = 4;

% place uncleanable object in scene
%-------------------------------------------------------------------------- 
d    = 1 & (xyz(:,1) > 2) & (xyz(:,1) < 8);
d    = d & (xyz(:,2) > (Ny - 10)) & (xyz(:,2) < (Ny - 2));
d    = d & (xyz(:,3) < 6);
W(d) = 7;

% place walls in scene
%-------------------------------------------------------------------------- 
W(1, :,:) = 6;                             % wall
W(Nx,:,:) = 6;                             % wall
W(:,1 ,:) = 6;                             % wall
W(:,Ny,:) = 6;                             % wall

% place supports in scene
%--------------------------------------------------------------------------
c = find(W(:,:,1) ~= 6);
c = c(randperm(numel(c),4));
for i = 1:numel(c)
    [ind] = spm_index([Nx,Ny],c(i));
    d     = 1 & (xyz(:,1) == ind(1));
    d     = d & (xyz(:,2) == ind(2));
    d     = d & (xyz(:,3) < Nz);
    W(d)  = 6;
end

% Now specify the hidden states at each location (a MAP)
%--------------------------------------------------------------------------
for i = 1:Nx
    for j = 1:Ny
        for k = 1:Nz
            D{i,j,k} = sparse(W(i,j,k),1,1,Nc,1);
        end
    end
end

return

function spm_behaviour(MDP)
% display posteriors
% FORMAT spm_show_x(x,Nx,Ny)
% a{g} - likelihood tensors
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% sensor's location and observations
%--------------------------------------------------------------------------
global PG

Nr    = PG.Nr*PG.Nr;                       % number of rays
N     = prod([PG.Nx,PG.Ny,PG.Nz]);         % number of locations
Nu    = 4;                                 % number of domain factors
X     = MDP.s(1:Nu,:);                     % location
L     = (1:N) + Nu;                        % state indices

% location of hidden factors in rectangular coordinates
%--------------------------------------------------------------------------
XYZ   = spm_combinations([PG.Nx,PG.Ny,PG.Nz]);


% illustrate behaviour
%--------------------------------------------------------------------------
for t = 1:MDP.T

    % true scene
    %----------------------------------------------------------------------
    Q     = {};
    for f = 1:N
        Q{f,1} = sparse(MDP.s(f + Nu,t),1,1,PG.Nc,1);
    end
    subplot(2,2,2), hold off
    spm_show_x(Q,PG.Nx,PG.Ny,PG.Nz), hold on
    title('True scene')

    % inferred scene
    %----------------------------------------------------------------------
    Q     = MDP.X(L);
    for f = 1:N
        Q{f} = Q{f}(:,t);
    end
    subplot(2,2,1), hold off
    spm_show_x(Q,PG.Nx,PG.Ny,PG.Nz), hold on
    title(sprintf('Inferred scene (t = %i)',t))

    % where the agent thinks it is
    %----------------------------------------------------------------------
    for f = 1:Nu
        [~,j] = max(MDP.X{f}(:,t));
        Y(f)  = j;
    end
    plot(Y(2),Y(1),    '*y','MarkerSize',32,'LineWidth',4)

    % true location and orientation of agent
    %----------------------------------------------------------------------
    plot(X(2,t),X(1,t),'*w','MarkerSize',32,'LineWidth',4)
    
    % Has the agent crashed
    %----------------------------------------------------------------------
    f  = find(ismember(XYZ,X(1:3,t)','rows'));
    if MDP.s(f + Nu,t) > 1
        plot(X(2,t),X(1,t),'*r','MarkerSize',64,'LineWidth',4)
    end

    % where the agent has been
    %----------------------------------------------------------------------
    plot(X(2,1:t),X(1,1:t),'.w','MarkerSize',16)

    % where the agent looks
    %----------------------------------------------------------------------
    j     = MDP.j(:,t);
    g     = MDP.id.ge;
    for i = g(g <= Nr)
        for k = 1:numel(j{i})
            x = spm_index([PG.Nx,PG.Ny,PG.Nx],j{i}(k) - Nu);
            plot(x(2),x(1),'ow','MarkerSize',x(3))
        end
    end

    % 3D view: true scene
    %======================================================================
    subplot(4,2,6), hold off
    for s = 2:PG.Nc
        i = find(MDP.s(L,t) == s);
        c = spm_colour(sparse(s,1,1,PG.Nc,1));
        plot3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'.','MarkerSize',32,'Color',c)
        hold on
    end
    
    % agent location
    %----------------------------------------------------------------------
    plot3(X(1,t),X(2,t),X(3,t),'*w','MarkerSize',16,'LineWidth',2), hold on
    plot3(X(1,t),X(2,t),2     ,'.w','MarkerSize',48,'LineWidth',2), hold on
    plot3([1,1]*X(1,t),[1,1]*X(2,t),[2 PG.Nz],':w', 'LineWidth',2), hold on
    axis image, a  = axis;

    % 3D view :inferred
    %======================================================================
    for f = 1:N
        [m,i] = max(Q{f});
        if m > 2/PG.Nc
            q(f) = i;
        else
            q(f) = 1;
        end
    end
    subplot(4,2,5), hold off
    for s = 2:PG.Nc
        i = find(q == s);
        c = spm_colour(sparse(s,1,1,PG.Nc,1));
        plot3(XYZ(i,1),XYZ(i,2),XYZ(i,3),'.','MarkerSize',32,'Color',c)
        hold on
    end
    
    % agent location
    %----------------------------------------------------------------------
    plot3(X(1,t),X(2,t),X(3,t),'*w','MarkerSize',16,'LineWidth',2), hold on
    plot3(X(1,t),X(2,t),2     ,'.w','MarkerSize',48,'LineWidth',2), hold on
    plot3([1,1]*X(1,t),[1,1]*X(2,t),[2 PG.Nz],':w', 'LineWidth',2), hold on
    axis image, axis(a)


    % what the agent sees (depth)
    %----------------------------------------------------------------------
    g     = 1:Nr;
    depth = MDP.o(g,t);
    depth = reshape(depth,PG.Nr,PG.Nr);
    subplot(4,3,10), hold off
    imagesc(1 - depth), axis image
    title('Depth')
    dmin  = min(depth(:));

    % what the agent sees (state or semantic class)
    %----------------------------------------------------------------------
    g     = (1:Nr) + Nr;
    state = MDP.O(g,t);
    state = reshape(state,PG.Nr,PG.Nr);

    % find class
    %--------------------------------------------------------------
    I     = zeros(PG.Nr,PG.Nr,3);
    J     = zeros(PG.Nr,PG.Nr,3);
    for i = 1:PG.Nr
        for j = 1:PG.Nr
            I(i,j,:) = spm_colour(state{i,j});
            J(i,j,:) = I(i,j,:)/((depth(i,j) - dmin + 1));
        end
    end
    subplot(4,3,12), hold off
    imagesc(I), axis image
    title('Class')

    subplot(4,3,11), hold off
    imagesc(J), axis image
    title('Vision')

    % save movie
    %----------------------------------------------------------------------
    drawnow
    MOV(t) = getframe(gcf);

end

% Place movie in graphic subject
%--------------------------------------------------------------------------
set(gcf,'Userdata',{MOV,2})
set(gcf,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

return




