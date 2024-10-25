function ADEM_physics
% % This routine illustrates the application of Variational Laplace to the
% triple-estimation (DEM) problem; namely inferring the latent states and
% causes of a dynamical system while learning the parameters and
% precisions. The example is chosen to highlight complementary ways of
% dealing with random (unknown) variables that change slowly. Here, these
% variables include the attribute of various objects (which are inferred as
% latent causes) and the parameters underwriting their dynamics (which are
% learned). The example in question is a simple Atari-like game but
% expressed in continuous (generalised) coordinates of motion.
% 
% To illustrate how one can equip DEM with action (i.e., implement an
% elementary form of active inference) it considers the problem of using a
% (red) bat to repel a (blue) ball, in the presence of other (green)
% objects. The dynamics underlying this behaviour are based upon
% straightforward classical mechanics, in which interactions amongst
% objects are parameterised with coupling parameters that have to be
% learned.
% 
% In what follows, we first set up the generative process and accompanying
% model. By exposing the scheme to the initial conditions (that evolve very
% slowly) the location and attributes of objects are then inferred to
% provide precise posteriors over latent causes of (visual) input. These
% posteriors and then used to learn the (coupling) parameters of the
% generative model that capture the interactions between the bat and ball.
% These dynamics are learned in the absence of action. Finally, we enable
% action such that the inference and learning scheme controls the ball by
% fulfilling its predictions under what it has previously learned. This
% results in the emergence of a nonequilibrium steady-state, here,
% quasi-periodic orbits that inherit from the physics of the generative
% model
% 
% The likelihood mapping from latent states generates visual input using a
% simple convolution model, in which each (Gaussian) object of an inferred
% size is convolved with (Gaussian) receptive fields. For biomimetic
% reasons, these receptive fields are larger at the periphery of vision.
% And can be regarded as an implementation of spatial basis functions
% mapping from object location to the (TrueColour) amplitude of place-coded
% responses: c.f., retinotopic map or a downsampled image after
% (anisotropic) convolution with receptive field kernels. The contributions
% of different wavelengths are generated following a softmax operator that
% stands in for a normalisation implicit in Weber's law. Despite these
% nonlinearities, model inversion can still proceed using local
% linearisation as in extended Kalman-Bucy filters. It should be noted,
% however, that changing any of the model parameters will preclude
% convergence, illustrating the failure Variational Laplace under brittle,
% nonlinear models. On the other hand, with the right parameterisation, the
% ensuing scheme illustrates the simplicity of specifying and inverting
% such models with just one inversion scheme (here, spm_ADEM).
%__________________________________________________________________________
 
% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

% Generative process and model
%==========================================================================
% First, set up the generative process and model by specifying the number
% of objects and sensory (visual) RGB channels
%--------------------------------------------------------------------------
rng(2);                                     % to ensure convergence
global RGB

% sensors
%--------------------------------------------------------------------------
Np     = 4;                                 % number of sensors
No     = 3;                                 % number of objects
RGB.Np = Np;                                % number of receptive fields
RGB.ij = spm_combinations([Np,Np]);         % their indices
RGB.xy = 8*(RGB.ij - mean(1:Np))/Np;        % position and
RGB.s  = sum(RGB.xy.^2,2)/4;                % size (variance)

RGB.action = 0;                             % switch off active inference

% parameters
%==========================================================================

% first object is repelled by second object
%--------------------------------------------------------------------------
P.k      = zeros(No,No);                    % coupling parameters
P.k(1,2) =  8;                              % 1st object attracted by 2nd
P.k(2,1) = -1;                              % 2nd object repelled by 1st

% parameters (of potential) for control object
%--------------------------------------------------------------------------
P.x  = [-3, 0];                             % equilibrium point
P.u  = [64, 0;                              % curvature
        0, 1];
       
% and parameters of a potential function of position for remaining objects
%--------------------------------------------------------------------------
P.sx  = zeros(1,2);                         % equilibrium point
P.su  = eye(2,2);                           % curvature

% specify priors (pE,pC) over parameters
%--------------------------------------------------------------------------
pE    = P;                                  % prior expectation
pC    = spm_zeros(P);                       % prior covariances
pE.k  = spm_zeros(P.k);                     % prior expectation (coupling)
pC.k  = (1 - eye(No,No))*exp(16);           % prior covariances (coupling)

% states (dynamics and features: i.e., what and where)
%--------------------------------------------------------------------------
for i = 1:No

    % object-specific states
    %----------------------------------------------------------------------
    x(i).x = randn(1,2);                    % position   
    x(i).u = randn(1,2);                    % velocity

    % object-specific attributes
    %----------------------------------------------------------------------
    v(i).c = randn(1,3);                    % log TrueColor
    v(i).s = randn;                         % log size (variance)

end

% first object
%--------------------------------------------------------------------------
for i = 1

    % object-specific states
    %----------------------------------------------------------------------
    x(i).x = P.x;
    x(i).u = [0, 0];

    % object-specific attributes
    %----------------------------------------------------------------------
    v(i).c = [-3, 0, 0];
    v(i).s = 0;

end

% priors over initial states and causes
%--------------------------------------------------------------------------
px    = spm_zeros(x);
pv    = spm_zeros(v);
for i = 1:No

    % object-specific states: initalise near centre of scene
    %----------------------------------------------------------------------
    px(i).x = randn(1,2)/32;             

end

% generative process
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = @spm_fx_aphysics;                   % dynamics
G(1).g  = @spm_gx_aphysics;                   % observation model
G(1).pE = P;                                  % parameters
G(1).x  = x;                                  % initial states
G(1).V  = exp(16);                            % precision (observation)
G(1).W  = exp(16);                            % precision (state noise)
G(1).U  = exp(4);                             % precision (action)

% second level
%--------------------------------------------------------------------------
G(2).v  = v;                                  % initial causes
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);                            % precision (latent causes)

% generative model
%==========================================================================
M(1).E.s  = 1;                               % smoothness
M(1).E.n  = 4;                               % order of 
M(1).E.d  = 0;                               % generalised motion
M(1).E.dt = 1/128;                           % time bin
M(1).E.nE = 1;                               % number of E-steps

% level 1: dynamics (where)
%--------------------------------------------------------------------------
M(1).f  = @spm_fx_physics;                   % physics
M(1).g  = @spm_gx_physics;                   % prediction
M(1).x  = px;                                % initial states
M(1).pE = pE;                                % parameters
M(1).V  = exp(8);                            % error precision (noise)
M(1).W  = exp(4);                            % error precision (state)

% level 2: what
%--------------------------------------------------------------------------
M(2).v  = pv;                                % inital causes
M(2).V  = exp(0);                            % with imprecise priors
 
% generate and invert - inference (D-Step)
%==========================================================================
% Using this model, we can now present the initial conditions (very slowly)
% to allow reasonably precise inferences about the latent causes; namely,
% the size and colour of objects based upon their inferred locations. This
% can be regarded as disentangling an (almost) static scene to find the
% best explanation in terms of a small number of objects with (almost)
% invariant characteristics.
%--------------------------------------------------------------------------
Nt      = 64;                                % length of data sequence
DEM.G   = G;
DEM.M   = M;
DEM.C   = repmat(spm_vec( v),[1 Nt]);        % latent causes
DEM.U   = repmat(spm_vec(pv),[1 Nt]);        % priors over causes

% state estimate (inference)
%--------------------------------------------------------------------------
DEM     = spm_ADEM(DEM);
spm_DEM_qU(DEM.qU,DEM.pU)

% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Segmentation'); clf
spm_dem_physics_movie(DEM)

% generate and invert - learning (E and M-steps)
%==========================================================================
% We can now use the posterior over the latent causes as a prior over the
% attribute of objects. This enables us to present a short episode of
% dynamics (i.e., gameplay) to learn the underlying parameters; here,
% restricted to the coupling parameters
%--------------------------------------------------------------------------

% Bayesian belief updating (latent causes)
%--------------------------------------------------------------------------
px  = DEM.qU.x{1}(:,end);
pv  = DEM.qU.v{2}(:,end);
px  = spm_unvec(px,x);
pv  = spm_unvec(pv,v);

M(1).E.dt = 1/2;                             % time bin
M(1).E.nE = 8;                               % number of E-steps

M(1).x  = px;                                % intial states
M(1).pE = pE;                                % prior expectation 
M(1).pC = pC;                                % and covariance (parameters)
M(1).gE = 0;                                 % prior (log-precision)
M(1).V  = exp(8);                            % error precision (noise)
M(1).W  = exp(0);                            % error precision (state)

% level 2: what
%--------------------------------------------------------------------------
M(2).v  = pv;                                % inital causes
M(2).V  = exp(8);                            % with precise priors
 
% Place process, model, and priors over latent causes in DEM structure
%--------------------------------------------------------------------------
Nt      = 32;                                % length of data sequence
DEM.G   = G;
DEM.M   = M;
DEM.C   = repmat(spm_vec( v),[1 Nt]);
DEM.U   = repmat(spm_vec(pv),[1 Nt]);

% Inference and learning
%--------------------------------------------------------------------------
DEM     = spm_ADEM(DEM);

% state estimates (inference)
%--------------------------------------------------------------------------
spm_figure('GetWin','States'); clf
spm_DEM_qU(DEM.qU,DEM.pU)

% parameter estimates (learning)
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameters'); clf
spm_DEM_qP(DEM.qP,DEM.pP)

% precision estimates (learning)
%--------------------------------------------------------------------------
spm_figure('GetWin','Precisions'); clf
spm_DEM_qH(DEM.qH)

% Graphics (latent and predicted outcomes)
%--------------------------------------------------------------------------
spm_figure('GetWin','Physics'); clf
spm_dem_physics_movie(DEM)


% Enable action to illustrate active inference
%==========================================================================
% Finally, we can now use the posterior expectations of the parameters to
% generate predictions in output space that action can fulfil by,
% effectively, minimising prediction errors with (merely reflexive) active
% inference
%--------------------------------------------------------------------------
RGB.action  = 1;

% update parameters and remove uncertainty
%--------------------------------------------------------------------------
DEM.M(1).pE = DEM.qP.P{1};                       % parameters
DEM.M(1).pC = spm_zeros(pC);                     % and covariance
DEM.M(1).gE = DEM.qH.g{1};                       % log precision
DEM.M(1).gC = spm_zeros(1);                      % and covariance

Nt      = 512;                                   % length of data sequence
DEM.C   = repmat(spm_vec( v),[1 Nt]);
DEM.U   = repmat(spm_vec(pv),[1 Nt]);

% active inference
%--------------------------------------------------------------------------
DEM  = spm_ADEM(DEM);

% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Active Inference'); clf
spm_dem_physics_movie(DEM)

return


function [f] = spm_fx_physics(x,v,P)
% returns the flow (f) for an n-body problem
% FORMAT [f] = spm_fx_physics(x,v,P)
%
% x    - latent states
% v    - latent causes
% P    - parameters
%
% This subroutine is a specification of the dynamics; in terms of the
% motion (flow) of position and velocity. The flow of velocity this
% acceleration induced by classical forces. These forces depend upon a
% shared potential function of position (V) and interactions among massive
% objects (K) scaled by coupling parameters. The requisite functions are
% the derivative of the respective potentials.
%__________________________________________________________________________

V  = @(x) x;                             % potential gradient
K  = @(x) x.*exp(-x.^2/2);               % coupling function

% dynamics
%--------------------------------------------------------------------------
f     = x;
No    = numel(x);
for i = 1:No

    % effective mass
    %----------------------------------------------------------------------
    m      = 8*exp(v(i).s);

    % forces exerted by other objects
    %----------------------------------------------------------------------
    F     = 0;
    for j = 1:No

        if j == i

            % friction
            %--------------------------------------------------------------
            F = F - P.k(i,i)*x(i).u;

        else

            % interactions
            %--------------------------------------------------------------
            F = F + P.k(i,j)*K(x(j).x - x(i).x)/m;
        end
    end

    % velocity
    %----------------------------------------------------------------------
    f(i).x = x(i).u;

    % acceleration
    %----------------------------------------------------------------------
    if i > 1
        f(i).u = V(P.sx - x(i).x)*P.su/m + F;
    else
        f(i).u = V(P.x - x(i).x)*P.u/m + F - x(i).u/2;
    end

end

return

function [g] = spm_gx_physics(x,v,P)
% returns the obervations (g) for an n-body problem
% FORMAT [g] = spm_gx_physics(x,v,P)
%
% x    - latent states
% v    - latent causes
% P    - parameters
%
% This subroutine specifies the nonlinear likelihood mapping between the
% latent variables and observations (g). The outcomes are generated by
% first, involving a Gaussian function of position with a Gaussian
% convolution kernel. This kernel is specific to each sensor location and
% can be regarded as a receptive field (as specified in the global variable
% RGB). The (log) colours are then normalised such that the (Np x Np x 3)
% output tensor contains contributions of each colour from each object. For
% convenience (when showing the observations in image format) we work with
% absence of colours; such that white denotes no nearby object.
%__________________________________________________________________________
global RGB

% observation
%--------------------------------------------------------------------------
No    = numel(x);
g.x   = ones(RGB.Np,RGB.Np,3);
for i = 1:No

    % shape (convolution with receptive fields)
    %----------------------------------------------------------------------
    X   = minus(RGB.xy,x(i).x);
    s   = exp(v(i).s);
    s   = s^2 + RGB.s;
    O   = reshape(exp(-sum(X.^2,2)./(2*s)),RGB.Np,RGB.Np);

    % colour (normalization over colour channels)
    %----------------------------------------------------------------------
    col = spm_softmax(v(i).c'*4);
    for j = 1:3
        g.x(:,:,j)  = g.x(:,:,j).*(1 - O*col(j));
    end

end

return

function [f] = spm_fx_aphysics(x,v,a,P)
% returns the flow (f) for an n-body problem
% FORMAT [f] = spm_fx_physics(x,v,a,P)
%
% x    - latent states
% v    - latent causes
% a    - action (conrol) variables 
% P    - parameters
%
% This subroutine is called by the generative process and has the same
% functional form as the dynamics of the generative model, with the
% exception that it expects an additional input argument (a). If action is
% enabled, the forces on the first (controllable) object are specified by
% action, which is evaluated automatically as the action that minimises
% variational free energy (i.e., minimises prediction errors in observation
% space)
%__________________________________________________________________________
global RGB

% dynamics
%--------------------------------------------------------------------------
f  = spm_fx_physics(x,v,P);

if RGB.action

    % acceleration
    %----------------------------------------------------------------------
    f(1).u = [a(1), a(2)];

end

return

function [g] = spm_gx_aphysics(x,v,a,P)
% returns the obervations (g) for an n-body problem
% FORMAT [g] = spm_gx_physics(x,v,a,P)
%
% x    - latent states
% v    - latent causes
% a    - action (conrol) variables 
% P    - parameters
%__________________________________________________________________________

% observation
%--------------------------------------------------------------------------
g  = spm_gx_physics(x,v,P);

return

function spm_dem_physics_movie(DEM)
% creates a movie of n-body problem
% FORMAT spm_dem_physics_movie(DEM)
%
% DEM - DEM structure from simulation
%__________________________________________________________________________

% true states
%==========================================================================
subplot(3,2,1)
M = spm_dem_physics_plot(DEM.pU,DEM.G);

% set ButtonDownFcn
%--------------------------------------------------------------------------
h = gca;
set(h(1),'Userdata',{M,16})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

% observations
%==========================================================================
subplot(3,2,3)
M = spm_dem_physics_image(DEM.pU,DEM.G);

% set ButtonDownFcn
%--------------------------------------------------------------------------
h = gca;
set(h(1),'Userdata',{M,16})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

% movie inferred states
%==========================================================================
subplot(3,2,2)
M = spm_dem_physics_plot(DEM.qU,DEM.M);


% set ButtonDownFcn
%--------------------------------------------------------------------------
h = gca;
set(h(1),'Userdata',{M,16})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')

% predictions
%==========================================================================
subplot(3,2,4)
M = spm_dem_physics_image(DEM.qU,DEM.M);

% set ButtonDownFcn
%--------------------------------------------------------------------------
h = gca;
set(h(1),'Userdata',{M,16})
set(h(1),'ButtonDownFcn','spm_DEM_ButtonDownFcn')



function M = spm_dem_physics_plot(U,G)
% subroutine for the movie: for plotting states and causes
% FORMAT M = spm_dem_physics_plot(U,G)
%
% U - states and causes (in standard format)
% G - Generative process or model (in standard format)
%__________________________________________________________________________
global RGB

% true states
%==========================================================================
for t = 1:size(U.v{2},2)

    cla
    hold on

    % true states
    %----------------------------------------------------------------------
    x    = spm_unvec(U.x{1}(:,t),G(1).x);
    v    = spm_unvec(U.v{2}(:,t),G(2).v);

    % plot
    %----------------------------------------------------------------------
    for i = 1:numel(x)
        siz = 4 + exp(v(i).s);
        col = spm_softmax(-v(i).c'*4)';
        plot(x(i).x(2),x(i).x(1),'.','MarkerSize',siz*8,'color',col)
        axis([-1 1 -1 1]*8), axis square
    end

    % pixels: i.e., location of receptive field centres
    %----------------------------------------------------------------------
    plot(RGB.xy(:,1),RGB.xy(:,2),'.k','MarkerSize',2)
    
    % save
    %----------------------------------------------------------------------
    drawnow
    M(t) = getframe(gca);

end

function M = spm_dem_physics_image(U,G)
% subroutine for the movie: for creating coloured images of observations
% FORMAT M = spm_dem_physics_image(U,G)
%
% U - states and causes (in standard format)
% G - Generative process or model (in standard format)
%__________________________________________________________________________
% observations
%==========================================================================
for t = 1:size(U.v{1},2)

    % true states
    %--------------------------------------------------------------------------
    g  = spm_unvec(U.v{1}(:,t),G(1).v);

    % image
    %----------------------------------------------------------------------
    imagesc(g.x), axis square, axis xy
    drawnow

    % save
    %----------------------------------------------------------------------
    M(t) = getframe(gca);

end
