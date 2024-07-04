function MDP = DEM_drone_V
% Demo of domain factors in the setting of active vision
%__________________________________________________________________________
%
% This routine illustrates a kind of SLAM (simultaneous localistion and
% mapping) powered by active inference. The generative model is lightweight
% in terms of memory; thereby eluding the von Neumann bottleneck. It
% introduces a number of innovations. First, large tensors are replaced by
% function handles that return requisite indices and probability
% distributions on the fly (as a function of domain factors and states,
% respectively). This precludes sum-product operators in the generative
% process and takes the pressure off precomputed indices conditioned upon
% domain factors. Note that outcomes are specified probabilistically;
% enabling suitable uncertainty about sensor inputs (e.g., semantic
% classification of locations that are out of range).
%
% Second, action is now handled more realistically by selecting the action
% that maximises the accuracy under posterior predictions. In other
% words, action is explicitly determined based on the consequences in a
% subspace of (telemetry) outcomes.
% 
% Third, this necessitates the provision of telemetry in terms of position,
% pose and (lidar) depth along rays. Finally, the depth of each ray (i.e.,
% line of sight) is treated as a latent state that is directly informed by
% (lidar) measurements — and used as a domain factor to select the
% appropriate likelihood mapping. This means that, in the generative model,
% there is no need to model occlusion; this modelling is offloaded to the
% generative process using the functional form of the requisite likelihood
% mapping.
% 
% The above scheme allows for more rays of greater depth than previously
% afforded by tensor parameterisation. However, the price paid is compute
% time, which scales with the number of rays. Biomimetically, this can be
% regarded as an exhaustive saccadic search over all rays. This follows
% because the belief updating accumulates evidence from the likelihood of
% each ray sequentially. For simplicity, the likelihood mappings are
% assumed to be very precise; given the course graining afforded by
% discrete state space models.
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

Nx  = 32;                                % size of environment
Ny  = 32;                                % size of environment
Nz  = 8;                                 % size of environment
Nd  = 16;                                % depth of rays
Nc  = 8;                                 % number of classes
Nr  = 6;                                 % number of rays
Na  = 12;                                % number of drone angles

% global parameters
%--------------------------------------------------------------------------
global PG %#ok<GVMIS> 
PG.Nx  = Nx;
PG.Ny  = Ny;
PG.Nz  = Nz;
PG.Na  = Na;
PG.Nd  = Nd;
PG.Nr  = Nr;
PG.Nc  = Nc;

PG.phi = linspace(0,2*pi*(1 - 1/Na),Na); % drone angles
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
title('True environment'),  drawnow

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


% Likelihood mapping: given states along the line of sight (i.e., rays)
%==========================================================================
% Given the state of locations at various depths along a line of
% sight, specify the outcome in terms of depth and class; namely, the
% state of the first non-empty location, where the first state corresponds
% to empty space.
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
GA{Ng + Ng + 4} = eye(Na,Na);              % a - drone angles

% controlled outcomes (telemetry)
%--------------------------------------------------------------------------
ID.control = (1:4) + Ng + Ng;

% for every combination of (Nu) domain factors (i.e., id.ff) specify
% projective geometry
%--------------------------------------------------------------------------
Nu    = 4;
ID.ff = 1:Nu;
ID.fg = @spm_fg;

% repeat for 'model'
%==========================================================================


% Likelihood mapping: given states along the line of sight (i.e., rays)
%==========================================================================
% Given the state of locations at various depths along a line of
% sight, specifying the outcome in terms of depth and class; namely, the
% state of the first non-empty location, where the first state corresponds
% to empty space.
%--------------------------------------------------------------------------

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
A{Ng + Ng + 4} = eye(Na,Na);              % a - drone angles

% central and peripheral (visual) modalities for planning
%--------------------------------------------------------------------------
xyz  = spm_combinations([Nr,Nr]);
d    = (Nr - 1).^2/2;
cen  = find(sum(minus(xyz,[(Nr + 1)/2,(Nr + 1)/2]).^2,2) <  2);
per  = find(sum(minus(xyz,[(Nr + 1)/2,(Nr + 1)/2]).^2,2) >= d);

% uncomment to restrict modalities subtending EFE
%--------------------------------------------------------------------------
id.ge = unique([cen; (Ng + cen); per; (Ng + per)]');

% uncomment to remove depth modality
%--------------------------------------------------------------------------
% id.g = {Ng:(Ng + Ng + 4)};

% convert likelihood arrays to logical
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
% a - drone angle
% s - state at location 1
% ...
% s - state at location Nx*Ny*Nz = Nf
%--------------------------------------------------------------------------
Nu    = 4;                                 % # controllable factors
id.ff = 1:Nu;                              % domain factors (x,y,x,a)
id.fg = @spm_fg;                           % parents of A{g}


% Now specify transition priors (i.e., dynamics)
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
    end
    for j = 1:Nx
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
    end
    for j = 1:Ny
        if j > (Ny - d + 1)
            BY(:,j,i) = full(sparse((Ny - d + 1),1,1,Ny,1));
        end
    end

    % vertical and rotational moves
    %----------------------------------------------------------------------
    BZ(:,:,i) = full(spm_speye(Nz,Nz,u(i),2));
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
n     = [1/16, 1/16, 1/64, 0]/128;               % state noise
for f = 1:numel(n)
    b           = GB{f};
    c           = n(f);
    GB{f}(:,:,1) = b(:,:,1) + c*b(:,:,2) + c*b(:,:,1)^2;
    GB{f}(:,:,2) = b(:,:,2) + c*b(:,:,1) + c*b(:,:,3);
    GB{f}(:,:,3) = b(:,:,3) + c*b(:,:,2) + c*b(:,:,3)^2;
    GB{f}        = spm_dir_norm(GB{f});
end

% add depth factors to model
%==========================================================================

% priors over initial states
%--------------------------------------------------------------------------
B     = GB;                                      % state transitions
for f = 1:numel(B)
    D{f} = ones(size(B{f},1),1);
end

V    = fix([Nx/2,Ny/2,Nz,1]);                    % initial telemetry
D{1} = full(sparse(V(1),1,1,Nx,1));              % location (x)
D{2} = full(sparse(V(2),1,1,Ny,1));              % location (y)
D{3} = full(sparse(V(3),1,1,Nz,1));              % location (z)
D{4} = full(sparse(V(4),1,1,Na,1));              % drone angles


%% priors: (cost) C
%==========================================================================
% no modality contraints at this stage
%--------------------------------------------------------------------------
for g = 1:numel(A)
    C{g} = spm_dir_norm(ones(size(A{g},1),1));
end

% In addition, specify constraints in latent sub-space (cif)
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
mdp.T  = T;                       % numer of moves
mdp.A  = A;                       % likelihood probabilities
mdp.B  = B;                       % transition probabilities
mdp.C  = C;                       % prior constraints
mdp.D  = D;                       % initial priors
mdp.N  = N;                       % policy depth
mdp.U  = U;                       % controllable actions

mdp.GA = GA;                      % likelihood probabilities
mdp.GB = GB;                      % transition probabilities
mdp.GU = U;                       % controllable actions

mdp.ID = ID;                      % process domains
mdp.id = id;                      % model domains
mdp.s  = s;                       % inital states

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
spm_behaviour(MDP)

return


%--------------------------------------------------------------------------
fprintf('How many 3-things did you find?\n\n')
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


% Solve - an example of search and rescue
%==========================================================================
% In the second task, the agent has to find something (i.e., class 3) and
% maintain surveillance over it. To simulate this, we specify a preference
% class 3 in the class or attribute modalities
%--------------------------------------------------------------------------
c     = spm_softmax(sparse(3,1,4,Nc,1));
for g = (1:Ng) + Ng
    C{g} = c;
end

MDP   = mdp;
MDP.C = C;
MDP.T = 16;
MDP   = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Search and rescue');
spm_behaviour(MDP)

%--------------------------------------------------------------------------
fprintf('I''ve found a 3-thing\n\n')
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
MDP.D((1:Nf) + Nu) = MAP;

% Specify intial (home) and final hidden (goal) states (hid)
%--------------------------------------------------------------------------

% find empty intial (s0) and final (sT) states
%--------------------------------------------------------------------------
for s = 1:prod([Nx,Ny,Nz])
    XYZ(s,:) = spm_index([Nx,Ny,Nz],s);
end
s          = XYZ(W == 1,:);
[~,d]      = min(sum(abs(minus(s,[Nx - 4,Ny - 4,1])),2)); sT = s(d,:);
[~,d]      = min(sum(abs(minus(s,[4,4,1])),2));           s0 = s(d,:);


MDP.s(1:4) = [s0, 1];
MDP.D{1}   = full(sparse(s0(1),1,1,Nx,1));   % inital state
MDP.D{2}   = full(sparse(s0(2),1,1,Ny,1));   % inital state
MDP.D{3}   = full(sparse(s0(3),1,1,Nz,1));   % inital state
MDP.D{4}   = full(sparse(    1,1,1,Na,1));   % inital state

MDP.id.hid = sT(:);                          % final  state
MDP        = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Wayfinding'); clf
spm_behaviour(MDP)
subplot(2,2,1), hold on, plot(sT(1),sT(2),'ow','MarkerSize',48)

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
MDP.D((1:Nf) + Nu) = MAP;

% active inference
%--------------------------------------------------------------------------
MDP.T = 8;
MDP   = spm_MDP_VB_XXX(MDP);

spm_figure('GetWin','Orientation');
spm_MDP_VB_trial(MDP);

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

% and preclude if not empty
%--------------------------------------------------------------------------
D   = Q(1,:) > (1 - 1/32);
D   = reshape(full(D),PG.Nx,PG.Ny,PG.Nz);

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
% s(4) - drone angles
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
% the drone can occupy. For example, if the drone can point in six
% directions, the field of view subtends an angle of 60°. Furthermore, we
% assume this field of view points downwards such that the drone can what
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
% GA{Ng + Ng + 4} = eye(Na,Na);              % a - drone angles
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

% get centre and angles of drone
%--------------------------------------------------------------------------
origin = [s(1:3)];
angle  = PG.phi(s(4));

% get sampled locations for this modality
%--------------------------------------------------------------------------
theta  = PG.FOV;                            % ray angles (60 degree FOV)
d      = spm_index([PG.Nr,PG.Nr],g);              % ray (line of sight) angle
ph     = angle - theta(d(2));               % azimuthal angle
th     = theta(d(1)) - theta(1);            % polar angle

% for each depth along ray (i.e., line of sight)
%--------------------------------------------------------------------------
for d = 1:PG.Nd

    % find nearest location in the latent state space
    %----------------------------------------------------------------------
    %     xyz    = spm_spherical2rectang(origin,ph,th,d);
    %     x      = uint8(max(min(xyz(1),Nx),1));
    %     y      = uint8(max(min(xyz(2),Ny),1));
    %     z      = uint8(max(min(xyz(3),Nz),1));
    %     j      = spm_sub2ind([Nx,Ny,Nz],x,y,z);
    %     par(d) = j + Nu;

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
MAP = [0 0 0;                     % (1) black (void)
    000 204 204;                  % (2) blue  (water)
    000 153 000;                  % (3) green (grass)
    051 102 000;                  % (4) green (foliage)
    204 204 000;                  % (5) green (pillar)
    255 229 204;                  % (6) pink  (person)
    204 102 000;                  % (7) brick (brick)
    255 000 000]/255;             % (8) red   (target)
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

function [D,W] = spm_drone_unity(Nx,Ny,Nz,Nc)
% creates an outdoor scene based on size of the environment

% Create a random environment with several objects of different classes
%--------------------------------------------------------------------------
% (1) black (void)
% (2) blue  (water)
% (3) green (grass)
% (4) green (foliage)
% (5) brown (pillar)
% (6) pink  (person)
% (7) brick (brick)
% (8) red   (target)

xyz      = spm_combinations([Nx,Ny,Nz]);
W        = ones(Nx,Ny,Nz);                % empty space
W(:,:,1) = 2;                             % water base

% place foliage in scene
%--------------------------------------------------------------------------
d    = sqrt(sum(minus(xyz(:,[1 2]),[Nx/2,Ny/3]).^2,2)) > Nx/4;
d    = d & (xyz(:,3) < 3);
W(d) = 3;                             % land (grass)

% place foliage in scene
%-------------------------------------------------------------------------- 
c = find(W(:,:,2) == 3);
c = c(randperm(numel(c),10));
for i = 1:numel(c)
    [ind] = spm_index([Nx,Ny],c(i));
    d     = sum(minus(xyz,[ind(1),ind(2),4 + rand]).^2,2);
    d     = sqrt(d) < Nx/(9 + rand);
    d     = logical(d.*(rand(size(d)) > .6));
    W(d)  = 4;
end

% place building in scene
%-------------------------------------------------------------------------- 
d    = 1 & (xyz(:,1) > 2) & (xyz(:,1) < 6);
d    = d & (xyz(:,2) > (Ny - 10)) & (xyz(:,2) < (Ny - 2));
d    = d & (xyz(:,3) < 6);
W(d) = 7;

% place pillars in scene
%--------------------------------------------------------------------------
c = find(W(:,:,2) == 3);
c = c(randperm(numel(c),4));
for i = 1:numel(c)
    [ind] = spm_index([Nx,Ny],c(i));
    d     = 1 & (xyz(:,1) == ind(1));
    d     = d & (xyz(:,2) == ind(2));
    d     = d & (xyz(:,3) < (Nz - rand));
    W(d)  = 5;
end

% place people in scene
%-------------------------------------------------------------------------- 
d    = 1 & (xyz(:,1) > 5) & (xyz(:,1) < 8);
d    = d & (xyz(:,2) > (Ny - 10)) & (xyz(:,2) < (Ny - 2));
d    = logical(d.*(rand(size(d)) > .5));
d    = d & (xyz(:,3) == 3);
W(d) = 6;

% place target in pond
%--------------------------------------------------------------------------
c     = find(W(:,:,1) == 2);
c     = c(randperm(numel(c),1));
[ind] = spm_index([Nx,Ny],c);
d     = 1 & (xyz(:,1) == ind(1));
d     = d & (xyz(:,2) == ind(2));
d     = d & (xyz(:,3) == 2);
W(d)  = 8;

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

% drone's location and observations
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
    spm_show_x(Q',PG.Nx,PG.Ny,PG.Nz), hold on
    title(sprintf('Inferred scene (t = %i)',t))

    % where the drone thinks it is
    %----------------------------------------------------------------------
    for f = 1:Nu
        [~,j] = max(MDP.X{f}(:,t));
        Xt(f) = j;
    end
    plot(Xt(2), Xt(1), '*y','MarkerSize',32,'LineWidth',4)

    % true location and orientation of drone
    %----------------------------------------------------------------------
    plot(X(2,t),X(1,t),'*w','MarkerSize',32,'LineWidth',4)
    
    % Has the drone crashed
    %----------------------------------------------------------------------
    f  = find(ismember(XYZ,X(1:3,t)','rows'));
    if MDP.s(f + Nu,t) > 1
        plot(X(2,t),X(1,t),'*r','MarkerSize',64,'LineWidth',4)
    end

    % where the drone has been
    %----------------------------------------------------------------------
    plot(X(2,1:t),X(1,1:t),'.w','MarkerSize',16)

    % where the drone looks
    %----------------------------------------------------------------------
    j     = MDP.j(:,t);
    for i = 1:Nr
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
    
    % drone location
    %----------------------------------------------------------------------
    plot3(X(1,t),X(2,t),X(3,t),'*w','MarkerSize',16,'LineWidth',2), hold on
    axis image, a  = axis;

    % 3D view :inferred
    %======================================================================
    for f = 1:N
        [m,i] = max(Q{f});
        if m > 1 - 1/16
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
    
    % drone location
    %----------------------------------------------------------------------
    plot3(X(1,t),X(2,t),X(3,t),'*k','MarkerSize',16,'LineWidth',2), hold on
    axis image, axis(a)


    % what the drone sees (depth)
    %----------------------------------------------------------------------
    g     = 1:Nr;
    depth = MDP.o(g,t);
    depth = reshape(depth,PG.Nr,PG.Nr);
    subplot(4,3,10), hold off
    imagesc(1 - depth), axis image
    title('Depth')
    dmin  = min(depth(:));

    % what the drone sees (state or semantic class)
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




