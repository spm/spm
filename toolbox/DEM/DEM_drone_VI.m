function MDP = DEM_drone_VI
% Demo of domain factors in renormalising models
%__________________________________________________________________________
%
% This demonstration routine builds upon DEM_drone_VI by finessing the
% generative model and, crucially, using fast structure learning to a
% symbol a scale free (renormalising) generative model: i.e., an RGM. In
% brief, this furnishes a generative model of a known scene that includes
% events. For example, people moving around the scene. In more detail, an
% RGM has learned automatically from latent states generated from a process
% model with no drone movement. The ensuing RGM then generates sequences of
% state transitions at the lowest level. To complete the RGM, the
% generative model encoding the lowest level likelihood mapping between
% latent states and (drone position and pose-dependent) outcomes replaces
% the first level of the RGM, with suitable adjustments to match the
% outputs of the second level RGM to the children in the first level (in a
% way that the paths of latent states of the first level are generated).
%
% Effectively, this equips the generative model with empirical priors that
% inherit from loan dynamics in the environment; enabling the agent to
% track moving things to resolve uncertainty about where things (e.g.,
% people) will appear next. Note that this demonstration effectively places
% knowledge about the scene at hand in the parameters of the likelihood
% mapping from the second to the first level. In other words, this
% demonstration is concerned with navigating in a known environment, where
% knowledge is installed both about things that do not move and dynamics.
% From the perspective of the latent states encoding the state of each
% location, this means that some locations never change, while an empty
% location can change to be occupied by the class of things that can move,
% and vice versa. The temporal or serial dependencies on these transitions
% encode movement and object permanence as one ascends the hierarchical
% generative model that features Scale transformations in both space and
% time.
%
% In this example, the empirical priors afforded by the deeper levels of
% the RGM applied to and only two the latent states associated to each
% location. This leaves the state of the drone (that would be informed by
% telemetry) as controllable states to illustrate planning as inference at
% the lowest level (using explicit action in the usual way). This means
% that the planning is only constrained or informed by posterior
% predictions about the next state of the world.
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

Nx  = 32;                                % size of environment
Ny  = 32;                                % size of environment
Nz  = 8;                                 % size of environment
Nd  = 16;                                % depth of rays
Nc  = 8;                                 % number of classes
Nr  = 16;                                % number of rays
Na  = 12;                                % number of drone angles
No  = 1;                                 % number of animate objects
Nu  = 4;                                 % # of control states
Ns  = 6;                                 % numer of FOV

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
PG.Nu  = Nu;
PG.No  = No;

PG.phi = linspace(0,2*pi*(1 - 1/Na),Na); % drone angles
PG.FOV = linspace(-pi/Ns,pi/Ns,Nr);        % ray angles (60 degree FOV)


% Get environment (inanimate scene)
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
Np    = size(XYZ,1);                       % # of places
Ng    = Nr*Nr;                             % # of rays
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
ID.control =  Ng + Ng + (1:Nu);

% for every combination of (Nu) domain factors (i.e., id.ff) specify
% projective geometry in functional form (see: spm_fg)
%--------------------------------------------------------------------------
ID.ff = 1:PG.Nu;
for o = 1:PG.No
    ff    = Nu + Np + (o - 1)*4 + (1:3);
    ID.ff = [ID.ff ff];
end
ID.fg = @spm_fg;

% repeat for 'model'
%==========================================================================

% depth and vision modalities: model
%--------------------------------------------------------------------------
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
xyz    = spm_combinations([Nr,Nr]);
d      = (Nr - 1).^2/2;
PG.cen  = find(sum(minus(xyz,[(Nr + 1)/2,(Nr + 1)/2]).^2,2) <  2);
PG.per  = find(sum(minus(xyz,[(Nr + 1)/2,(Nr + 1)/2]).^2,2) >= d);

% uncomment to restrict modalities subtending EFE
%--------------------------------------------------------------------------
id.ge = unique([PG.cen; (Ng + PG.cen); PG.per; (Ng + PG.per)]');

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

% latent states:
%==========================================================================
% x - drone location (x)
% y - drone location (y)
% z - drone location (z)
% a - drone angle
% x - agent location (x)
% y - agent location (y)
% z - agent location (z)
% s - state of agent (c)
% ...
% s - state at location 1 (c)
% ...
% s - state at location Nx*Ny*Nz = Nl
%--------------------------------------------------------------------------
id.ff = 1:PG.Nu;                           % domain factors (x,y,x,a)
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

    % floor contraints on z
    %----------------------------------------------------------------------
    BZ(:,:,i) = full(spm_speye(Nz,Nz,u(i),2));
    for j = 1:Nz
        if j < 3
            BZ(:,j,i) = full(sparse(3,1,1,Nz,1));
        end
    end

    % rotational moves
    %----------------------------------------------------------------------
    BA(:,:,i) = full(spm_speye(Na,Na,u(i),1));

end

% control factors
%==========================================================================
GB    = {BX,BY,BZ,BA};                           % control factors

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



% append latent place factors (process)
%--------------------------------------------------------------------------
Bp    = logical(eye(Nc,Nc));
for f = 1:Np                             
    GB{end + 1} = Bp;
end

% state transitions for model: these will be replaced later when mdp by
% empirical priors furnished by an RGM
%--------------------------------------------------------------------------
B      = GB; 

% append process transitions for animate object (e.g., person)
%--------------------------------------------------------------------------
u      = 1;
d      = 1;
e      = 12;
Bx(:,:,u)   = full(spm_speye(Nx,Nx,-1,1));
Bx(:,1:d)   = 0;   Bx(d + 1,1:d) = 1;
Bx(:,e:end) = 0; Bx(d + 1,e:end) = 1;

By(:,:,u) = full(spm_speye(Ny,Ny));
Bz(:,:,u) = full(spm_speye(Nz,Nz));
Bc(:,:,u) = full(spm_speye(Nc,Nc));
for  o = 1:PG.No
    GB = [GB,{Bx},{By},{Bz},{Bc}];
end

% priors over initial states
%==========================================================================
for f = 1:numel(B)
    D{f} = ones(size(B{f},1),1);
end

V    = fix([Nx/2,Ny/2,Nz,2*Na/3]);               % initial telemetry
D{1} = full(sparse(V(1),1,1,Nx,1));              % location (x)
D{2} = full(sparse(V(2),1,1,Ny,1));              % location (y)
D{3} = full(sparse(V(3),1,1,Nz,1));              % location (z)
D{4} = full(sparse(V(4),1,1,Na,1));              % drone angles


% priors: (cost) C
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
GU    = zeros(1,numel(GB));
GU(u) = 1;

% true (process) states
%--------------------------------------------------------------------------
P     = fix([4,Ny/3,3,6]);        % initial object states
P     = repmat(P,[1,PG.No]);
s     = [V(:); W(:); P(:)];

% MDP Structure, specifying T epochs of active vision
%==========================================================================
mdp.A  = A;                       % likelihood probabilities
mdp.B  = B;                       % transition probabilities
mdp.C  = C;                       % prior constraints
mdp.D  = D;                       % initial priors
mdp.N  = 0;                       % policy depth
mdp.U  = U;                       % controllable actions
mdp.id = id;                      % model domains

mdp.GA = GA;                      % likelihood probabilities
mdp.GB = GB;                      % transition probabilities
mdp.GU = GU;                      % controllable actions
mdp.ID = ID;                      % process domains
mdp.s  = s;                       % inital states

% true (training) states
%--------------------------------------------------------------------------
P     = fix([1,Ny/3,3,6]);        % initial object states
P     = repmat(P,[1,PG.No]);
s     = [V(:); W(:); P(:)];

% generative process
%--------------------------------------------------------------------------
GDP.T  = 64;                      % time-steps
GDP.A  = GA;                      % likelihood probabilities
GDP.B  = GB;                      % transition probabilities
GDP.U  = GU*0;                    % controllable actions
GDP.id = ID;                      % process domains
GDP.s  = s;                       % inital states

% generate a dynamic scene
%--------------------------------------------------------------------------
PDP   = spm_MDP_generate(GDP);
Q     = cell(Np,GDP.T);
for t = 1:GDP.T

    % background
    %----------------------------------------------------------------------
    for f = 1:Np
        Q{f,t} = sparse(PDP.s(f + Nu,t),1,1,PG.Nc,1);
    end

    % add objects
    %----------------------------------------------------------------------
    for o = 1:PG.No
        s = PDP.s(Nu + Np + (o - 1)*4 + (1:4),t);
        f = spm_sub2ind([PG.Nx,PG.Ny,PG.Nz],s(1),s(2),s(3));
        Q{f,t} = sparse(s(4),1,1,PG.Nc,1);
    end
end

subplot(2,2,2), hold off
spm_show_x(Q,PG.Nx,PG.Ny,PG.Nz)
title('True scene')


% RG structure learning
%==========================================================================
tic,
MDP = spm_faster_structure_learning(Q,[Nx,Ny,Nz],[1 64]);
t   = toc;

spm_figure('GetWin',sprintf('Paramters: level %i',numel(MDP))); clf
spm_MDP_params(MDP{end}), drawnow

fprintf('Structure learning time: %i secs.\n',round(t))

% normalise Dirichlet coutns and remove a and b fields
%--------------------------------------------------------------------------
MDP = spm_mdp_a2A(MDP);

% replace first level likelihood mapping (and augment with telemetry)
%--------------------------------------------------------------------------
mdp.id.D = cell(1,numel(mdp.B));
mdp.id.E = cell(1,numel(mdp.B));
for f = 1:numel(MDP{1}.B)

    % expand B
    %----------------------------------------------------------------------
    for u = 1:size(MDP{1}.B{f},3)
        b = MDP{1}.A{f}*MDP{1}.B{f}(:,:,u)*MDP{1}.A{f}';
        mdp.B{Nu + f}(:,:,u) = spm_dir_norm(b);
    end

    % and parents
    %----------------------------------------------------------------------
    g                = MDP{1}.id.D{f};
    a                = MDP{1}.A{f}*MDP{2}.A{g};
    MDP{2}.A{g}      = spm_dir_norm(a);
    mdp.id.D(Nu + f) = MDP{1}.id.D(f);
    mdp.id.E(Nu + f) = MDP{1}.id.E(f);

end

% Create deep recursive model
%--------------------------------------------------------------------------
RDP   = [{mdp} MDP(2:end)];
RDP   = spm_mdp2rdp(RDP);

% Solve - an example of surveillance
%==========================================================================
% In this instance, we expect to see the drone orientate to dynamic parts
% of the scene and keep an eye on moving objects (e.g., people). The
% generative process underwriting the RGM means that the drone a priori
% believes there are repeating events, where each event there is a simple
% orbit (e.g., someone entering the scene and disappearing in the water, to
% re-emerge the top of the scene).
%--------------------------------------------------------------------------
RDP.T    = 16;
tic; RDP = spm_MDP_VB_XXX(RDP); t = toc;

% Illustrate recursive model
%--------------------------------------------------------------------------
spm_figure('GetWin','RGM'); clf
spm_show_RGM(RDP);

% illustrate scene construction and exploration
%--------------------------------------------------------------------------
spm_figure('GetWin','Active inference');
QDP    = spm_Q2mdp(RDP.Q,1);
QDP.id = mdp.id;
fprintf('Compute time: %i ms/update\n',round(1000*t/size(QDP.O,2)))

spm_behaviour(QDP)


% Solve - an example of wayfinding with priors
%==========================================================================
% In this example, we Illustrate pathfinding using inductive inference
% under precise prior beliefs about latent states (i.e., cognitve map).
%
% This simulation is identical to the previous illustration; however, now
% the agent is equipped with delivery or goal states, which it navigates
% to, while maintaining its situational awareness (i.e., with and eye
% dynamic scene components).
%--------------------------------------------------------------------------
pdp        = mdp;
pdp.D((1:Np) + Nu) = MAP;

% Specify intial (home) and final hidden (goal) states (hid)
%--------------------------------------------------------------------------
sT         = [28,12, 3];                     % final  state
pdp.id.hid = sT(:);                          % final  state

RDP        = [{pdp} MDP(2:end)];
RDP        = spm_mdp2rdp(RDP);
RDP.T      = 4;
RDP        = spm_MDP_VB_XXX(RDP);

spm_figure('GetWin','Wayfinding'); clf
QDP    = spm_Q2mdp(RDP.Q,1);
QDP.id = pdp.id;

spm_behaviour(QDP)
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
d   = find(s > 1,1,'first');             % first nonempty class
if numel(d)
    O(d,1)         = 1;                  % depth of first occluder
else
    O(PG.Nd + 1,1) = 1;                  % nothing in range
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
    O = p.*cumprod([1 q(1:end - 1)]);      % P(first occupied)
    O = O(:);                              % P(outcome)
    O(PG.Nd + 1,1) = max(1 - sum(O),0);
    O = O/sum(O);

elseif nargin == 1

    % generate from states
    %----------------------------------------------------------------------
    O = spm_depth(P);

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
    d       = A_depth(P);                 % depth of first occluder
    q       = spm_cat(P);                 % P(state)
    q(1,:)  = 0;                          % P(state | occupied)
    q(1,PG.Nd + 1) = 1;                   % out of range

    q = spm_dir_norm(q);                  % normalise
    O = q*d;                              % P(outcome)

elseif nargin == 1

    % generate from states
    %----------------------------------------------------------------------
    O = spm_state(P);

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
Np  = PG.Nx*PG.Ny*PG.Nz;
Q   = spm_cat(Q((1:Np) + PG.Nu));

% and preclude if not empty
%--------------------------------------------------------------------------
D   = Q(1,:) > (1 - 1/32);
D   = reshape(full(D),PG.Nx,PG.Ny,PG.Nz);

% (location) factors
%--------------------------------------------------------------------------
hif = 1:3;

return

function par = spm_fg(g,s)
% projective geometry
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

Nu    = PG.Nu;                               % domain factors
Ng    = PG.Nr*PG.Nr;                         % number of rays
Np    = PG.Nx*PG.Ny*PG.Nz;                   % number of places

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
d      = spm_index([PG.Nr,PG.Nr],g);        % ray (line of sight) angle
ph     = angle - theta(d(2));               % azimuthal angle
th     = theta(d(1)) - theta(1);            % polar angle

% location of (animate) object
%--------------------------------------------------------------------------
No    = (numel(s) - Nu)/3;
jo    = 0;
for o = 1:No
    so    = s(Nu + (o - 1)*3 + (1:3));           % place
    fo(o) = Nu + Np + o*4;                       % class
    jo(o) = spm_sub2ind([PG.Nx,PG.Ny,PG.Nz],so(1),so(2),so(3));
end

% for each depth along ray (i.e., line of sight)
%--------------------------------------------------------------------------
par   = zeros(1,0);
for d = 1:PG.Nd

    % find nearest location in the latent state space
    %----------------------------------------------------------------------
    xyz = ceil(spm_spherical2rectang(origin,ph,th,d));
    x   = uint8(xyz(1));
    y   = uint8(xyz(2));
    z   = uint8(xyz(3));
    if x > 0 && x <= PG.Nx && y > 0 && y <= PG.Ny && z > 0 && z <= PG.Nz

        % within scene location
        %------------------------------------------------------------------
        j      = spm_sub2ind([PG.Nx,PG.Ny,PG.Nz],x,y,z);

        % is this location occupied by an (animate) object?
        %------------------------------------------------------------------
        par(d) = j + Nu;
        for o = 1:No
            if j == jo(o)
                par(d) = fo(o);
            end
        end
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
% x  - states
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% illustrate images
%--------------------------------------------------------------------------
u     = 1 - 1/2;
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
    drawnow

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
% FORMAT spm_behaviour(MDP)
% MDP.s  - states
% MDP.o  - outcomes
% MDP.X  - states   (posterior)
% MDP.O  - outcomes (posterior)
% MDP.j  - parents of outcomes
% MDP.id - edges
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% drone's location and observations
%--------------------------------------------------------------------------
global PG

Ng    = PG.Nr*PG.Nr;                       % number of rays
Np    = prod([PG.Nx,PG.Ny,PG.Nz]);         % number of places
Nu    = PG.Nu;                             % number of domain factors
X     = MDP.s(1:Nu,:);                     % drone location
L     = (1:Np) + Nu;                       % place indices

% location of hidden factors in rectangular coordinates
%--------------------------------------------------------------------------
XYZ   = spm_combinations([PG.Nx,PG.Ny,PG.Nz]);

% illustrate behaviour
%--------------------------------------------------------------------------
for t = 1:size(MDP.s,2)

    % true scene
    %----------------------------------------------------------------------
    Q     = {};
    for f = 1:Np
        Q{f,1} = sparse(MDP.s(f + Nu,t),1,1,PG.Nc,1);
    end

    % add objects
    %----------------------------------------------------------------------
    try
        for o = 1:PG.No
            s = MDP.s(Nu + Np + (o - 1)*4 + (1:4),t);
            f = spm_sub2ind([PG.Nx,PG.Ny,PG.Nz],s(1),s(2),s(3));
            Q{f,1} = sparse(s(4),1,1,PG.Nc,1);
        end
    end

    subplot(2,2,2), hold off
    spm_show_x(Q,PG.Nx,PG.Ny,PG.Nz), hold on
    title('True scene')

    % inferred scene
    %----------------------------------------------------------------------
    Q     = MDP.X(L);
    for f = 1:Np
        Q{f} = Q{f}(:,t);
    end
    subplot(2,2,1), hold off
    spm_show_x(Q,PG.Nx,PG.Ny,PG.Nz), hold on
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
    g     = MDP.id.ge;
    for i = g(g <= Ng)
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

    % add object
    %----------------------------------------------------------------------
    try
        for o = 1:PG.No
            s = MDP.s(Nu + Np + (o - 1)*4 + (1:4),t);
            f = spm_sub2ind([PG.Nx,PG.Ny,PG.Nz],s(1),s(2),s(3));
            c = spm_colour(sparse(s(4),1,1,PG.Nc,1));
            plot3(XYZ(f,1),XYZ(f,2),XYZ(f,3),'.','MarkerSize',32,'Color',c)
        end
    end

    % drone location
    %----------------------------------------------------------------------
    plot3(X(1,t),X(2,t),X(3,t),'*w','MarkerSize',16,'LineWidth',2), hold on
    axis image, a  = axis;

    % 3D view :inferred
    %======================================================================
    for f = 1:Np
        [m,i] = max(Q{f});
        if m > 1 - 1/8
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
    g     = 1:Ng;
    depth = MDP.o(g,t);
    depth = reshape(depth,PG.Nr,PG.Nr);
    subplot(4,3,10), hold off
    imagesc(1 - depth), axis image
    title('Depth')
    dmin  = min(depth(:));

    % what the drone sees (state or semantic class)
    %----------------------------------------------------------------------
    g     = (1:Ng) + Ng;
    state = MDP.O(g,t);
    state = reshape(state,PG.Nr,PG.Nr);

    % find class
    %----------------------------------------------------------------------
    I     = zeros(PG.Nr,PG.Nr,3);
    J     = zeros(PG.Nr,PG.Nr,3);
    for i = 1:PG.Nr
        for j = 1:PG.Nr
            I(i,j,:) = spm_colour(state{i,j});
            d        = (depth(i,j) - dmin)/PG.Nd;
            J(i,j,:) = I(i,j,:)/(4*d + 1);
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




