function ADEM_intercept_II
% Radar tracking and interception in generalised coordinates of motion
%__________________________________________________________________________
% This iteration of the demo considers learning about the missile as a
% prelude to online inference and tracking. The key move here is to observe
% the missile (for 100 seconds) and estimate the (time-invariant)
% parameters that govern its trajectory. These include things like the
% speed, source, target and any parameters encoding evasive manoeuvres.
% Crucially, these parameters are time-invariant and can be recovered using
% dynamic expectation maximisation (DEM) from an observed timeseries. After
% a sufficient period of observation, the parameters can then be used as
% hidden causes for active inference and tracking as implemented in ADEM.
% Note that these hidden causes can change from moment to moment because
% causes are time-dependent variables. In this example, an initial missile
% approach simulated until the missile is within 1 km of the line of
% defence. It then reverts to a linear trajectory towards the target and
% the gimbal infrared is engaged. 
%
% In addition, this demo allows for multiple radars and a re-factoring of
% the subroutines that subtend the priors and likelihood mappings for the
% generative process and model (where the behaviour of the missile is
% shared by an observation model for off-line parameter estimation and the
% generative model used for online tracking).
% 
% This demo accommodates conditional termination of the simulations
% specified in DEM.break.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

rng(1)

% hidden states, causes and parameters
%==========================================================================
% x    - hidden states:
%  x.m.x    - position (missile) (Km)
%  x.d.x    - position (drone)   (Km)
%
% v    - hidden causes:
%   v  - position (target), evasion coefficients, speed, etc
%
% g    - sensor modalities (for one radar): 
%   g(1) - bearing of missile (radar)   (radians)
%   g(2) - log range of missile (radar) (log(Km))
%   g(3) - position x (IMU)             (Km)
%   g(4) - position y (IMU)             (Km)
%   g(5) - bearing of missile (IR)      (radians)
%   g(6) - angular velocity   (IR)      (radians/s)
%--------------------------------------------------------------------------

% missile classes and models
%==========================================================================
% First, set up different classes of missiles in terms of their
% characteristics; i.e., maximum speed, source and target locations, nature
% of evasive motion, et cetera. We then specify a particular missile
% picking a combination of parameters for subsequent simulation.
%
% After this, we then specify a series of models and processes that will be
% used to classify the whistle in terms of its parameters (using DCM) and
% then track the missile by transferring the parameters to hidden causes of
% an online (ADEM) scheme.
%--------------------------------------------------------------------------
speed = num2cell([170 135 100 65]/60/60);     % max speed
weave = {[0,0],[1/2,0],[1/2,1/4]};            % evasive
point = num2cell([-2 -1/2 1]);                % target

% parameters
%--------------------------------------------------------------------------
P.r.x = [4, -4;-2, -1];                       % position  (radars)
P.d.s = 85/60/60;                             % max speed (drone)  per sec
P.m.s = speed{2};                             % max speed (missile) per sec
P.u   = [4; 12];                              % source location (Km)
P.v   = point{1};                             % target location (Km)
P.w   = weave{3};                             % evasion weights

% hidden states (initial)
%--------------------------------------------------------------------------
x.m.x = P.u;                                  % position (missile)
x.d.x = [0;-4];                               % position (drone)

% precisions of state noise (assumed by model)
%--------------------------------------------------------------------------
W.m.x = exp([16; 16]);                        % position (missile)
W.d.x = exp([16; 16]);                        % position (drone)

% precisions of observation noise
%--------------------------------------------------------------------------
n      = size(P.r.x,2);                       % number of radar
V.rad  = zeros(2,n) + exp(16);                % range & bearing (radar)
V.imu  = zeros(2,1) + exp(2);                 % position (IMU)
V.b    = 0;                                   % bearing (IR)
V.dbdt = 0;                                   % velocity (IR)


% precisions (for action): restricted to IMU signals
%--------------------------------------------------------------------------
U     = spm_zeros(V);                         % action precision
U.imu = U.imu + exp(4);                       % position (IMU)


% generative models
%==========================================================================

% order of generalised motion (and smoothness)
%--------------------------------------------------------------------------
dt     = 1;                                   % time step (seconds)
E.s    = 1/2;                                 % smoothness (seconds)
E.n    = 3;                                   % order (states)
E.d    = 0;                                   % order (causes)
E.dt   = dt;                                  % time step

% generative model for ADEM
%==========================================================================

% level 1: dynamics of hidden states
%--------------------------------------------------------------------------
M(1).f  = @fx_M;                              % fx (model)
M(1).g  = @gx_M;                              % gx (model)
M(1).x  = x;                                  % hidden states
M(1).V  = spm_vec(V);                         % error precision (signal)
M(1).W  = spm_vec(W);                         % error precision (dynamics)
M(1).E  = E;

% level 2: hidden cause (missile parameters)
%--------------------------------------------------------------------------
M(2).v  = spm_zeros(P);                       % hidden causes
M(2).V  = exp(32);                            % hidden causes

% generative process
%==========================================================================

% first level
%--------------------------------------------------------------------------
G(1).f  = @fx_G;                              % fx (process)
G(1).g  = @gx_G;                              % gx (process)
G(1).x  = x;                                  % hidden states
G(1).V  = exp(32);                            % error precision (errors)
G(1).W  = exp(32);                            % error precision (motion)
G(1).U  = spm_vec(U);                         % error precision (action)
G(1).pE = P;                                  % parameters

% second level
%--------------------------------------------------------------------------
G(2).v  = P;                                  % hidden causes
G(2).a  = [0; 0];                             % action
G(2).V  = exp(32);

% Observer model (and process)
%==========================================================================

% priors over parameters
%--------------------------------------------------------------------------
pE      = P;                                  % prior expectations
pE.m.s  = 100/60/60;                          % max speed (missile) per sec
pE.u    = pE.u;                               % source location (Km)
pE.v    = 0*pE.v;                             % target location (Km)
pE.w    = 0*pE.w;                             % evasion weights

pC      = spm_zeros(pE);                      % prior covariance
pC.m.s  = 1/256;                              % max speed (missile) per sec
pC.u    = pC.u + 0;                           % source location (Km)
pC.v    = pC.v + 64;                          % target location (Km)
pC.w    = pC.w + 1;                           % evasion weights

% first level
%--------------------------------------------------------------------------
O(1).f  = @fx_O;                              % fx (process)
O(1).g  = @gx_O;                              % gx (process)
O(1).x  = x.m.x;                              % hidden states
O(1).V  = exp(32);                            % error precision (errors)
O(1).W  = exp(32);                            % error precision (motion)
O(1).pE = pE;                                 % parameters
O(1).pC = pC;                                 % parameters
O(1).E  = E;

% second level
%--------------------------------------------------------------------------
O(2).v  = 0;

% Assemble DEM structure for missile observation
%==========================================================================
N       = 100/dt;                             % length of time series (s)
U       = zeros(1,N);                         % priors over latent causes

% generate (radar) observations under true parameters for 100 seconds
%--------------------------------------------------------------------------
DEM     = spm_DEM_generate(O,U,P);

% Estimate parameters using DEM
%--------------------------------------------------------------------------
DEM.M(1).E.nE = 16;

DEM.M(1).V = spm_vec(V.rad);
DEM.M(1).W = exp(8);
DEM        = spm_DEM(DEM);

% show posterior estimate of trajectory
%--------------------------------------------------------------------------
spm_figure('GetWin','Missile classification'); clf

subplot(3,1,1)
Q  = DEM.qP.P{1};
px = DEM.pU.x{1};
qx = DEM.qU.x{1};
plot(px(1,:),px(2,:),'.'),  hold on
plot(qx(1,:),qx(2,:),'.r')
plot([P.u(1),P.v(1)],[P.u(2),0],':b')
plot([P.u(1),P.v(1)],[P.u(2),0],'ob')
plot([Q.u(1),Q.v(1)],[Q.u(2),0],':r')
plot([Q.u(1),Q.v(1)],[Q.u(2),0],'or')
for i = 1:n
    plot(P.r.x(1,i),P.r.x(2,i),'^k','MarkerSize',16), hold on
end
axis([-8,8,-4,12]), xlabel('Km'), ylabel('Km')
title('Initial trajectories'), axis square

% and illustrate posteriors over parameters
%--------------------------------------------------------------------------
spm_figure('GetWin','Parameters'); clf
spm_DEM_qP(DEM.qP,DEM.pP)

% illustrate Bayesian model reduction based on posteriors
%==========================================================================
% This section illustrates the use of Bayesian model comparison (i.e. BMR)
% to show how a missile can be classified in terms of its characteristics
% by evaluating the variational free energy (i.e., elbow) under reduced
% priors afforded by the kinds of missiles that could be encountered. One
% could use the reduced posteriors as the basis of subsequent simulation;
% however, here, we just illustrate the marginal posterior distribution
% over different missile characteristics.
%--------------------------------------------------------------------------
spm_figure('GetWin','Missile classification');

% number of missile characteristics
%--------------------------------------------------------------------------
Ns = numel(speed);
Nw = numel(weave);
Nt = numel(point);

% priors and posteriors
%--------------------------------------------------------------------------
pE = DEM.M(1).pE;
pC = DEM.M(1).pC;
qE = DEM.qP.P{1};
qC = DEM.qP.C;

% Evaluate reduced free energy for all combinations of characteristics
%--------------------------------------------------------------------------
S     = spm_combinations([Ns,Nw,Nt]);
F     = zeros(size(S,1),1);
for i = 1:size(S,1)

    % reduced priors
    %----------------------------------------------------------------------
    rE     = pE;
    rC     = pC/256;
    rE.m.s = speed{S(i,1)};
    rE.v   = point{S(i,3)};
    rE.w   = weave{S(i,2)};
    
    % ELBO
    %----------------------------------------------------------------------
    F(i) = spm_log_evidence(qE,qC,pE,pC,rE,rC);

end
p   = spm_softmax(F(:));
p   = reshape(p,Ns,Nw,Nt);
m   = spm_marginal(p);

% show marginal probabilities over different missiles
%--------------------------------------------------------------------------
subplot(3,3,4), bar(m{1}), title('Speed'), axis square
subplot(3,3,5), bar(m{2}), title('Weave'), axis square
subplot(3,3,6), bar(m{3}), title('Target'),axis square

subplot(3,1,3), imagesc(squeeze(spm_sum(p,2))), axis square
xlabel('Speed'),ylabel('Target'),title('Joint probability')

% Initial interception
%==========================================================================
% In the next section, we now transfer the posterior over parameters to the
% hidden causes of a generative model. This allows the in causes to be
% refined on the basis of further observations. Crucially, we relax the
% prior precision over the target location in case the missile changes
% direction.
%--------------------------------------------------------------------------

% update states
%--------------------------------------------------------------------------
M(1).x.m.x = DEM.qU.x{1}(:,end);
G(1).x.m.x = DEM.pU.x{1}(:,end);

% update causes
%--------------------------------------------------------------------------
Q       = DEM.qP.P{1};                        % empirical priors
P       = P;                                  % true causes

% and transfer posterior over parameters to hidden causes (U)
%--------------------------------------------------------------------------
qV      = spm_inv(qC);
i       = spm_fieldindices(Q,'v');
qV(i,:) = 0;
qV(:,i) = 0;

% level 2: hidden causes
%--------------------------------------------------------------------------
M(2).v  = Q;                                  % initial causes
M(2).V  = qV;                                 % precision

% Assemble DEM structure
%--------------------------------------------------------------------------
clear DEM
N       = 400/dt;                             % length of time series (s)
DEM.G   = G;                                  % generative process
DEM.M   = M;                                  % generative model
DEM.C   = repmat(spm_vec(P),1,N);             % latent causes
DEM.U   = repmat(spm_vec(Q),1,N);             % priors over latent causes

% Solve or integrate
%--------------------------------------------------------------------------
DEM.break = 'pu.x{1}(2) < 1';
DEM       = spm_ADEM_UNITY(DEM);

% show results of Bayesian filtering
%--------------------------------------------------------------------------
spm_figure('GetWin','ADEM'); clf
spm_DEM_qU(DEM.qU)
subplot(2,2,1), title({'Sensor signals', 'prediction and error'},'FontSize',16)
subplot(2,2,2), title({'State estimation', 'hidden states'},'FontSize',16)

% create movie of radar tracking and intial approach
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_graphics_intercept(DEM)


% continue simulation at half the update speed (i.e., 250ms updates)
%==========================================================================
% Finally, update hidden states and causes and switch on the precision of
% the gimbal IR for close proximity interception.
%--------------------------------------------------------------------------

% update hidden and latent states
%--------------------------------------------------------------------------
DEM.M(1).x = spm_unvec(DEM.qU.x{1}(:,end),DEM.M(1).x);
DEM.M(1).v = spm_unvec(DEM.qU.v{1}(:,end),DEM.M(1).v);
DEM.M(2).v = spm_unvec(DEM.qU.v{2}(:,end),DEM.M(2).v);
DEM.G(1).x = spm_unvec(DEM.pU.x{1}(:,end),DEM.G(1).x);
DEM.G(2).v = spm_unvec(DEM.pU.v{2}(:,end),DEM.G(2).v);

% increase precision of IR signals
%--------------------------------------------------------------------------
V.b        = exp(0);
V.dbdt     = exp(0);                   
DEM.M(1).V = spm_vec(V);

% solve
%--------------------------------------------------------------------------
N     = 100/dt;                              % length of time series (s)
DEM.C = repmat(spm_vec(P),1,N);              % latent causes
DEM.U = repmat(spm_vec(Q),1,N);              % priors over latent causes

DEM.break = 'pu.x{1}(2) < 0';
FEM   = spm_ADEM_UNITY(DEM);                 % solve

% create movie of radar tracking and final approach
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_graphics_intercept(FEM,1)

return


% model functions (subroutines)
%__________________________________________________________________________

function f = fx_G(x,v,a,~)
% FORMAT f = fx_G(x,v,a,P)
% equations of motion for generative process
%--------------------------------------------------------------------------
%  x.m.x    - position (missile)
%  x.d.x    - position (drone)
%--------------------------------------------------------------------------
f     = x;

%  f.m.x - changes in position (missile)
%--------------------------------------------------------------------------
f.m.x = f_m(x.m.x,v);

%  f.d.x - changes in position (drone)
%--------------------------------------------------------------------------
f.d.x = a(1:2);

return


function f = fx_M(x,v,~)
% FORMAT f = fx_M(x,v,P)
% equations of motion for model
%--------------------------------------------------------------------------
%  x.m.x - position (missile)
%  x.d.x - position (drone)
%--------------------------------------------------------------------------
f     = x;

%  f.m.x - changes in position (missile)
%--------------------------------------------------------------------------
f.m.x = f_m(x.m.x,v);

%  f.d.x    - changes in position (drone)
%--------------------------------------------------------------------------
d     = [v.v; 0] - x.d.x;
r     = sqrt(d'*d);
d     = d/r;
f.d.x = erf(r*8)*v.d.s*d;

return

function f = fx_O(x,~,P)
% FORMAT f = fx_O(x,v,P)
% equations of motion for observation model
%--------------------------------------------------------------------------
%  x.m.x - position (missile)
%--------------------------------------------------------------------------

%  f.m.x - changes in position (missile)
%--------------------------------------------------------------------------
f = f_m(x,P);

return

function f = f_m(x,v)
% FORMAT f = f_m(x,v)
% equations of motion for model
%--------------------------------------------------------------------------
%  x  - position (missile)
%--------------------------------------------------------------------------

%  f  - changes in position (missile)
%--------------------------------------------------------------------------
v.v   = [v.v; 0];                               % target [x,y]
dx    = v.v - v.u;                              % origin to target
s     = sqrt(dx'*dx);                           % distance to origin
d     = dx/s;                                   % direction from origin
dm    = v.v - x;                                % missile to target
dv    = sqrt(dm'*dm);                           % distance to target
if dv > 1

    % add weaving (discrete cosine set)
    %----------------------------------------------------------------------
    r     = d\dm;                               % projection
    r     = r/s;                                % relative distance
    o     = [-d(2); d(1)];
    for i = 1:numel(v.w)
        d = d + v.w(i)*cos(pi*r*8*i)*o;
    end

else

    % final (linear) approach
    %----------------------------------------------------------------------
    d  = dm;

end


% velocity (if not over line of defence)
%--------------------------------------------------------------------------
if x(2) > 0
    f = v.m.s*d/sqrt(d'*d);
else
    f = 0*d;
end

return

function g = gx_G(x,v,~,P)
% FORMAT g = gx_G(x,v,a,P)
% observation function (process)
%--------------------------------------------------------------------------
g = gx_M(x,v,P);

return

function g = gx_M(x,v,P)
% FORMAT g = gx_M(x,v,P)
% observation function (model)
%--------------------------------------------------------------------------
%   g(1) - bearing of missile (radar)   (radians)
%   g(2) - log range of missile (radar) (log(Km))
%   g(3) - position x (IMU)             (Km)
%   g(4) - position y (IMU)             (Km)
%   g(5) - bearing of missile (IR)      (radians)
%   g(6) - angular velocity   (IR)      (radians/s)
%--------------------------------------------------------------------------
g          = gx_O(x.m.x,P,v);

%   g - position (IMU) 
%--------------------------------------------------------------------------
g(end + 1) = x.d.x(1);
g(end + 1) = x.d.x(2);

%   g - bearing of missile (IR gimble)
%--------------------------------------------------------------------------
d          = x.m.x - x.d.x;
b          = atan(d(1)/d(2));
g(end + 1) = b;

%   g - radial velocity (IR gimble)
%--------------------------------------------------------------------------
f          = f_m(x.m.x,v);
dbdt       = 1/(b^2 + 1)*(f(1)/d(2) - f(2)*d(1)/(d(2)^2));
g(end + 1) = dbdt;

return

function g = gx_O(x,~,P)
% FORMAT g = gx_O(x,v,P)
% observation function (model)
%--------------------------------------------------------------------------
%   g(1,i) - bearing of missile (radar)   (radians)
%   g(2,i) - log range of missile (radar) (log(Km))
%--------------------------------------------------------------------------
n     = size(P.r.x,2);
g     = zeros(2,n);
for i = 1:n

    %   g(i,1) - bearing of missile (radar)
    %--------------------------------------------------------------------------
    dx     = x - P.r.x(:,i);
    g(1,i) = atan(dx(1)/dx(2));

    %   g(i,2) - range of missile (radar)
    %--------------------------------------------------------------------------
    r      = sqrt(dx'*dx);
    g(2,i) = log(r);

end
g  = g(:);

return

% end of model functions
%__________________________________________________________________________


function spm_graphics_intercept(DEM,OPT)
% FORMAT spm_graphics_intercept(DEM,OPT)
%--------------------------------------------------------------------------
pP = DEM.pP.P{1};
px = DEM.pU.x{1};
pv = DEM.pU.v{2};
qx = DEM.qU.x{1};
qv = DEM.qU.v{2};
gv = DEM.qU.v{1};
T  = size(pv,2);
nr = size(pP.r.x,2);
if nargin < 2

    subplot(2,2,1)
    for t = 1:T

        % plot radar
        %------------------------------------------------------------------
        hold off,
        for i = 1:nr
            plot(pP.r.x(1,i),pP.r.x(2,i),'^k','MarkerSize',16), hold on
        end

        % plot missile, target and drone (actual)
        %------------------------------------------------------------------
        x = spm_unvec(px(:,t),DEM.G(1).x);
        v = DEM.G(1).pE.v;
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        plot(x.m.x(1),x.m.x(2),'.k','MarkerSize',32)
        plot(x.d.x(1),x.d.x(2),'*k','MarkerSize',16)
        plot(x.d.x(1),x.d.x(2),'ok','MarkerSize',8)
        plot(v(1),0,'+k','MarkerSize',24)

        x = spm_unvec(qx(:,t),DEM.M(1).x);
        v = spm_unvec(qv(:,t),DEM.M(2).v);

        % plot missile, target and drone (inferred)
        %------------------------------------------------------------------
        iv    = spm_fieldindices(DEM.M(2).v,'v');
        c     = DEM.qU.C{t}(iv,iv);
        c     = sqrt(c)*3;
        [i,j] = ellipsoid(v.v(1),0,1,c(1),1/4,0,32);

        fill(i(16,:)',j(16,:)',[1 1 1],'LineStyle',':','FaceColor','none','EdgeColor',[1 1 1]/8);
        plot(x.m.x(1),x.m.x(2),'.r','MarkerSize',16)
        plot(x.d.x(1),x.d.x(2),'*r','MarkerSize',8)
        plot(x.d.x(1),x.d.x(2),'or','MarkerSize',4)
        plot(v.v,0,'+r','MarkerSize',12)

        % plot gimble
        %------------------------------------------------------------------
        u(1) = x.d.x(1);
        u(2) = x.d.x(1) + sin(0);
        q(1) = x.d.x(2);
        q(2) = x.d.x(2) + cos(0);
        plot(u,q,'m')

        % plot bearing
        %------------------------------------------------------------------
        v   = -8;
        plot(g(nr*2 + 3)*8,v,'+m','MarkerSize',32)
        plot(g(nr*2 + 3)*8,v,'.m','MarkerSize',32)
        plot([-pi,pi],[v,v],':m')

        axis([-1 1 -1 1]*10), axis square
        xlabel('latitude (km)'), ylabel('longitude (km)'), title('Interception')
        

        % plot trajectories
        %------------------------------------------------------------------
        plot(px(1,1:t),px(2,1:t),':k')
        plot(px(3,1:t),px(4,1:t),':r')
        drawnow

        % save frame
        %------------------------------------------------------------------
        % M(t) = getframe(gca);

    end

    % save movie
    %----------------------------------------------------------------------
    % set(gca,'Userdata',[])
    % set(gca,'Userdata',{M,8})
    % set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

else

    subplot(2,2,2), hold off
    m     = zeros(2,T);
    d     = zeros(2,T);
    for t = 1:T

        % plot missile, target and drone (actual)
        %------------------------------------------------------------------
        x = spm_unvec(px(:,t),DEM.G(1).x);
        v = DEM.G(1).pE.v;
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        % distance between missile and drone
        %------------------------------------------------------------------
        m(:,t) = x.m.x;
        d(:,t) = x.d.x;

        hold off, plot(v(1),0,'+k','MarkerSize',24), hold on
        plot(x.m.x(1),x.m.x(2),'.k','MarkerSize',64)
        plot(x.d.x(1),x.d.x(2),'*k','MarkerSize',32)
        plot(x.d.x(1),x.d.x(2),'ok','MarkerSize',16)

        x = spm_unvec(qx(:,t),DEM.M(1).x);
        v = spm_unvec(qv(:,t),DEM.M(2).v);
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        % plot missile, target and drone (inferred)
        %------------------------------------------------------------------
        plot(x.m.x(1),x.m.x(2),'.r','MarkerSize',16)
        plot(x.d.x(1),x.d.x(2),'*r','MarkerSize',8)
        plot(x.d.x(1),x.d.x(2),'or','MarkerSize',4)
        plot(v.v,0,'+r','MarkerSize',12)

        % plot gimble
        %------------------------------------------------------------------
        u(1) = x.d.x(1);
        u(2) = x.d.x(1) + sin(0)/8;
        q(1) = x.d.x(2);
        q(2) = x.d.x(2) + cos(0)/8;
        plot(u,q,'m')

        axis([-1 1 -1 1] + [v.v v.v 0 0]), axis square
        xlabel('latitude (km)'), ylabel('longitude (km)'), title('Interception')

        % save frame
        %------------------------------------------------------------------
        drawnow
        M(t) = getframe(gca);

    end

    % save movie
    %----------------------------------------------------------------------
    set(gca,'Userdata',[])
    set(gca,'Userdata',{M,8})
    set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

    % plot trajectories
    %----------------------------------------------------------------------
    subplot(2,2,3), hold off
    m(1,:) = m(1,:) - v.v(1);
    m(2,:) = m(2,:) - 0;
    d(1,:) = d(1,:) - v.v(1);
    d(2,:) = d(2,:) - 0;

    plot(m(1,:)*1000,m(2,:)*1000,':k'), hold on
    plot(d(1,:)*1000,d(2,:)*1000,':r'), hold on
    axis([-1 1 -1 1]*10), axis square
    xlabel('latitude (m)'), ylabel('longitude (m)'), title('Trajectories (meters)')

    % plot distances (after interpolation)
    %----------------------------------------------------------------------
    im(1,:) = interp(m(1,:),1000);
    im(2,:) = interp(m(2,:),1000);
    id(1,:) = interp(d(1,:),1000);
    id(2,:) = interp(d(2,:),1000);

    dx = im - id;
    r  = sqrt(sum(dx.^2));

    subplot(2,2,4), hold off
    semilogy(r*1000,'.'), hold on
    semilogy(r*0 + 1e-0,'-.r')
    semilogy(r*0 + 1e-1,'-.r')
    str = sprintf('Distance (m) (min = %i cm)',fix(min(r)*100000));
    xlabel('time (ms)'), ylabel('meters'), title(str)
    axis square

end
