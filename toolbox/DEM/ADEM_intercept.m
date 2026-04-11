function ADEM_intercept
% Radar tracking and interception in generalised coordinates of motion
%__________________________________________________________________________
% This script illustrates the use of generalised filtering to solve the
% dual problem of radar tracking and interception. This is an elemental
% instance of (mere) active inference in which the priors over dynamics are
% realised through action. In this instance, the priors over latent states
% induce a dynamics with a particular fixed point; namely, a drone lies on
% the trajectory of a missile. It has an infrared sensor that points
% towards the oncoming missile and, crucially, the angular velocity of its
% gimbal mounted sensor is zero. The only situation in which these
% conditions are met is when the drone is on the path of an approaching
% missile. Writing down the equations of motion of the drone—with the above
% solution— effectively creates a generative model of a drone that, a
% priority, expects to intercept a missile. Conceptually, this can be read
% as an instance of the equilibrium point hypothesis for motor control in
% neurobiology. In economics, this would be regarded as a particular kind
% of prior; i.e., a heuristic. In control theory, this is a special
% instance of path integral control.
%
% In detail, the generative process (and model) involves three objects and
% a hidden cause. The first object is a missile travelling at a constant
% velocity (135 km/h towards a target or rendezvous point in some line of
% defence. This rendezvous point is the hidden cause. The second object is
% a drone that is attracted to the rendezvous point with a maximum speed
% (of 85 km/h). The third object is a (fixed location) radar that returns
% precise estimates of the missile’s range and bearing (in terms of log
% distance and radians, respectively). The drone is equipped with a
% gimbal-mounted IR camera that reports the bearing of the missile when
% sufficiently close. The drone’s IMU provides estimates of the drone’s
% velocity and current position, along with the angle and angular velocity
% of the gimbal. The active states of the drone effectively resolve IMU
% (c.f., proprioceptive) prediction errors based upon the predictions of
% the drone's dynamics and its gimbal.
%
% Intuitively, one can regard the hidden cause as setting the direction of
% travel for the missile and drone: to a common (rendezvous) point in a
% fixed line of defence. Initially, the radar tracking supplies a posterior
% estimate of the hidden cause and thereby attracts the drone to the
% intercept position. As the missile comes into the range of the drone’s
% infrared camera, posterior uncertainty about the hidden cause is
% resolved. In addition, the drone manoeuvres along the defence line to
% ensure its gimbal is pointing in the fixed direction of the oncoming
% missile.
%
% The inversion of the implicit continuous state-space model is in
% generalised coordinates of motion. To attain the requisite accuracy
% (i.e., centimetres) it is necessary to use relatively high orders of
% generalised motion (here, 8). This setup is very sensitive to the
% precisions of random fluctuations (and the relative estimates of these
% precisions used during generalised filtering or state estimation). In the
% example below, we have simplified the problem to 2 dimensions and have
% used very precise radar and IMU signals to demonstrate what could be
% achieved under realistic speeds and distances.
%
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

rng(0)

% hidden states, causes and parameters
%==========================================================================
% x    - hidden states:
%  x.m.x    - position (missile) (Km)
%  x.m.v    - velocity (missile) (Km/s)
%  x.d.x    - position (drone)   (Km)
%  x.d.v    - velocity (drone)   (Km/s)
%  x.g.x    - angle    (gimble)  (radians)
%  x.g.v    - velocity (gimble)  (radians/s)
%
% v    - hidden causes:
%   v  - position (target)
%
% g    - sensor modalities:
%   g(1) - bearing of missile (radar)   (radians)
%   g(2) - log range of missile (radar) (log(Km))
%   g(3) - velocity x (IMU)             (Km/s)
%   g(4) - velocity y (IMU)             (Km/s)
%   g(5) - bearing of missile (IR)      (radians)
%   g(6) - angle of gimble (IMU)        (radians)
%   g(7) - velocity gimble (IMU)        (Km/s)
%   g(8) - position x (IMU)             (Km)
%   g(9) - position y (IMU)             (Km)
%
% Parameters:
%   P.r.x - position (radar)
%   P.m.s - speed  (missile)
%   P.d.s - speed    (drone)
%--------------------------------------------------------------------------
P.r.x = [4;-2];                               % position  (radar)
P.m.s = 135/60/60;                            % max speed (missile)
P.d.s = 85/60/60;                             % max speed (drone)

% Initial hidden states
%--------------------------------------------------------------------------
x.m.x = [4; 8];                               % position (missile)
x.m.v = [0; 0];                               % velocity (missile)
x.d.x = [0;-4];                               % position (drone)
x.d.v = [0; 0];                               % velocity (drone)
x.g.x = 0;                                    % angle    (gimble)
x.g.v = 0;                                    % velocity (gimble)

% inital velocity (missile)
%--------------------------------------------------------------------------
dx    = [0; 0] - x.m.x;
x.m.v = P.m.s*dx/sqrt(dx'*dx);

% precisions of observation noise (assumed by model)
%--------------------------------------------------------------------------
V(1)  = exp(8);                               % bearing (radar)
V(2)  = exp(16);                              % range   (radar)
V(3)  = exp(8);                               % velocity x (IMU)
V(4)  = exp(8);                               % velocity y (IMU)
V(5)  = 0;                                    % bearing of missile (IR gimble)
V(6)  = exp(0);                               % angle of gimble
V(7)  = exp(0);                               % velocity (gimble)
V(8)  = exp(8);                               % position x (IMU)
V(9)  = exp(8);                               % position y (IMU)

% precisions of state noise (assumed by model)
%--------------------------------------------------------------------------
W.m.x = exp([16; 16]);                        % position (missile)
W.m.v = exp([16; 16]);                        % velocity (missile)
W.d.x = exp([16; 16]);                        % position (drone)
W.d.v = exp([16; 16]);                        % velocity (drone)
W.g.x = exp(16);                              % angle    (gimble)
W.g.v = exp(16);                              % velocity (gimble)


% generative model
%==========================================================================

% order of generalised motion (and smoothness)
%--------------------------------------------------------------------------
M(1).E.s  = 1;                                % smoothness
M(1).E.n  = 8;                                % order of
M(1).E.d  = 0;                                % generalised motion
M(1).E.dt = 1;                                % time step

% level 1: dynamics of hidden states
%--------------------------------------------------------------------------
M(1).f  = @fx_M;                              % fx (model)
M(1).g  = @gx_M;                              % gx (model)
M(1).x  = x;                                  % hidden states
M(1).V  = spm_vec(V);                         % error precision
M(1).W  = spm_vec(W);                         % error precision
M(1).pE = P;                                  % parameters


% level 2: hidden cause (rendezvous point)
%--------------------------------------------------------------------------
M(2).v  = [0; 0];                             % hidden causes
M(2).V  = [exp(-4); exp(16)];                 % [s.d.: about 8 km in x]

% generative process
%==========================================================================

% precisions (for action): restricted to IMU signals
%--------------------------------------------------------------------------
U(1)  = 0;                                    % bearing of missile (radar)
U(2)  = 0;                                    % range   of missile (radar)
U(3)  = exp(8);                               % velocity x (IMU)
U(4)  = exp(8);                               % velocity y (IMU)
U(5)  = 0;                                    % bearing of missile (IR gimble)
U(6)  = exp(8);                               % angle of gimble (IMU)
U(7)  = exp(8);                               % velocity gimble (IMU)
U(8)  = exp(8);                               % position x (IMU)
U(9)  = exp(8);                               % position y (IMU)

% precision of random fluctuations on sensor signals
%--------------------------------------------------------------------------
V(1)  = ((2*pi/360/10))^(-2);      % bearing (radar)   [sd: 1/10 degrees]
V(2)  = (1/1000/10)^(-2);          % range  (radar)    [sd: 1m at 10km]
V(3)  = (1/1000/10)^(-2);          % velocity x (IMU)  [sd: 1/10 m/s]
V(4)  = (1/1000/10)^(-2);          % velocity y (IMU)  [sd: 1/10 m/s]
V(5)  = ((2*pi/360/10))^(-2);      % bearing (IR)      [sd: 1/10 degrees]
V(6)  = ((2*pi/360/10))^(-2);      % angle of gimble   [sd: 1/10 degrees]
V(7)  = ((2*pi/360/10))^(-2);      % velocity (gimble) [sd: 1/10 degrees]
V(8)  = (1/1000/10)^(-2);          % position x (IMU)  [sd: 1/10 m]
V(9)  = (1/1000/10)^(-2);          % position y (IMU)  [sd: 1/10 m]

% precision of random fluctuations on latent states
%--------------------------------------------------------------------------
W.m.x = exp([28; 28]);                        % position (missile)
W.m.v = exp([16; 16]);                        % velocity (missile)
W.d.x = exp([28; 28]);                        % position (drone)
W.d.v = exp([16; 16]);                        % velocity (drone)
W.g.x = exp(28);                              % angle    (gimble)
W.g.v = exp(16);                              % velocity (gimble)

% first level
%--------------------------------------------------------------------------
G(1).f  = @fx_G;                              % fx (process)
G(1).g  = @gx_G;                              % gx (process)
G(1).x  = x;                                  % hidden states
G(1).V  = spm_vec(V);                         % error precision (errors)
G(1).W  = spm_vec(W);                         % error precision (motion)
G(1).U  = spm_vec(U);                         % error precision (action)
G(1).pE = P;                                  % parameters

% second level
%--------------------------------------------------------------------------
G(2).v  = [-2; 0];                            % hidden cause
G(2).a  = [0; 0; 0];                          % action
G(2).V  = exp(32);


% Assemble DEM structure
%==========================================================================
N       = 200;                                % length of time series (s)
DEM.G   = G;                                  % generative process
DEM.M   = M;                                  % generative model
DEM.C   = zeros(2,N);                         % latent causes
DEM.U   = zeros(2,N);                         % priors over latent causes
DEM.C(1,:) = -2;                              % target at x = -2, y = 0 Km

% Solve or integrate
%--------------------------------------------------------------------------
DEM     = spm_ADEM(DEM);

% show results of Bayesian filtering
%--------------------------------------------------------------------------
spm_figure('GetWin','ADEM'); clf
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(2,2,1), title({'Sensor signals', 'prediction and error'},'FontSize',16)
subplot(2,2,2), title({'State estimation', 'hidden states'},'FontSize',16)

% create movie of radar tracking and intial approach
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_graphics_intercept(DEM)


% continue simulation at half the update speed (i.e., 500ms updates)
%==========================================================================
DEM.M(1).E.dt = 1/2;                            % time step

% update hidden and latent states
%--------------------------------------------------------------------------
DEM.M(1).x = spm_unvec(DEM.qU.x{1}(:,end),DEM.M(1).x);
DEM.G(1).x = spm_unvec(DEM.pU.x{1}(:,end),DEM.G(1).x);
DEM.M(2).v = DEM.qU.v{2}(:,end);

% increase precision of IR signals and priors over hidden cause
%--------------------------------------------------------------------------
DEM.M(1).V(5,5) = exp(2);
DEM.M(2).V(1)   = exp(4);
DEM.U(1,:)      = DEM.qU.v{2}(1,end);

% solve
%--------------------------------------------------------------------------
DEM  = spm_ADEM(DEM);

% create movie of radar tracking and final approach
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_graphics_intercept(DEM,1)

return


function f = fx_G(x,v,a,P)
% FORMAT f = fx_G(x,v,a,P)
% equations of motion for generative process
%--------------------------------------------------------------------------
%  x.m.x    - position (missile)
%  x.m.v    - velocity (missile)
%  x.d.x    - position (drone)
%  x.d.v    - velocity (drone)
%  x.g.x    - angle    (gimble)
%  x.g.v    - velocity (gimble)
%--------------------------------------------------------------------------
f      = x;

%  f.m.x    - changes in position (missile)
%--------------------------------------------------------------------------
f.m.x  = x.m.v;

%  f.m.v    - changes in velocity (missile)
%--------------------------------------------------------------------------
dx     = v - x.m.x;
r      = sqrt(dx'*dx);
dv     = dx/r;
if dx(2) < 0
    f.m.v = P.m.s*dv - x.m.v;
else
    f.m.v = x.m.v - x.m.v;
end

%  f.d.x    - changes in position (drone)
%--------------------------------------------------------------------------
f.d.x  = x.d.v;

%  f.d.v    - changes in velocity (drone)
%--------------------------------------------------------------------------
f.d.v  = a(1:2);

%  f.g.x    - changes in angle    (gimble)
%--------------------------------------------------------------------------
f.g.x  = a(3);

%  f.g.v    - changes in velocity (gimble)
%--------------------------------------------------------------------------
f.g.v  = 0;

return


function f = fx_M(x,v,P)
% FORMAT f = fx_M(x,v,P)
% equations of motion for model
%--------------------------------------------------------------------------
%  x.m.x    - position (missile)
%  x.m.v    - velocity (missile)
%  x.d.x    - position (drone)
%  x.d.v    - velocity (drone)
%  x.g.x    - angle    (gimble)
%  x.g.v    - velocity (gimble)
%--------------------------------------------------------------------------
f     = x;

%  f.m.x    - changes in position (missile)
%--------------------------------------------------------------------------
dx    = v - x.m.x;
r     = sqrt(dx'*dx);
dv    = dx/r;
f.m.x = P.m.s*dv;

%  f.m.v    - changes in velocity (missile)
%--------------------------------------------------------------------------
f.m.v = x.m.v - x.m.v;

%  f.d.x    - changes in position (drone)
%--------------------------------------------------------------------------
f.d.x = x.d.v + [x.g.v; 0];

%  f.d.v    - changes in velocity (drone)
%--------------------------------------------------------------------------
dx    = v - x.d.x;
r     = sqrt(dx'*dx);
dv    = dx/r;
f.d.v = erf(r*8)*P.d.s*dv - x.d.v;

%  f.g.x    - changes in angle    (gimble)
%--------------------------------------------------------------------------
dx    = x.m.x - x.d.x;
b     = atan(dx(1)/dx(2));
f.g.x = b - x.g.x;

%  f.g.v    - changes in velocity (gimble)
%--------------------------------------------------------------------------
f.g.v = 0;

return


function g = gx_G(x,v,a,P)
% FORMAT g = gx_G(x,v,a,P)
% observation function (process)
%--------------------------------------------------------------------------
%   g(1) - bearing of missile (radar)
%   g(2) - range   of missile (radar)
%   g(3) - velocity x (IMU)
%   g(4) - velocity y (IMU)
%   g(5) - bearing of missile (IR)
%   g(6) - angle of gimble (IMU)
%   g(7) - velocity gimble (IMU)
%--------------------------------------------------------------------------
g      = zeros(9,1);

%   g(1) - bearing of missile (radar)
%--------------------------------------------------------------------------
dx     = x.m.x - P.r.x;
g(1)   = atan(dx(1)/dx(2));

%   g(2) - range of missile (radar)
%--------------------------------------------------------------------------
r      = sqrt(dx'*dx);
g(2)   = log(r);

%   g(3) - velocity (IMU)
%--------------------------------------------------------------------------
g(3:4) = x.d.v;

%   g(5) - bearing of missile (IR gimble)
%--------------------------------------------------------------------------
dx     = x.m.x - x.d.x;
g(5)   = atan(dx(1)/dx(2)) - x.g.x;

%   g(6) - angle of gimble (IMU)
%--------------------------------------------------------------------------
g(6)   = x.g.x;

%   g(7) - velocity gimble (IMU)
%--------------------------------------------------------------------------
g(7)   = x.g.v;

%   g(8:9) - (IMU)
%--------------------------------------------------------------------------
g(8:9) = x.d.x;


return

function g = gx_M(x,v,P)
% FORMAT g = gx_M(x,v,P)
% observation function (model)
%--------------------------------------------------------------------------
%   g(1) - bearing of missile (radar)
%   g(2) - range   of missile (radar)
%   g(3) - velocity x (IMU)
%   g(4) - velocity y (IMU)
%   g(5) - bearing of missile (IR)
%   g(6) - angle of gimble (IMU)
%   g(7) - velocity gimble (IMU)
%--------------------------------------------------------------------------
g      = zeros(9,1);

%   g(1) - bearing of missile (radar)
%--------------------------------------------------------------------------
dx     = x.m.x - P.r.x;
g(1)   = atan(dx(1)/dx(2));

%   g(2) - range of missile (radar)
%--------------------------------------------------------------------------
r      = sqrt(dx'*dx);
g(2)   = log(r);

%   g(3) - velocity (IMU)
%--------------------------------------------------------------------------
g(3:4) = x.d.v;

%   g(5) - bearing of missile (IR gimble)
%--------------------------------------------------------------------------
dx     = x.m.x - x.d.x;
g(5)   = atan(dx(1)/dx(2)) - x.g.x;

%   g(6) - angle of gimble (IMU)
%--------------------------------------------------------------------------
g(6)   = x.g.x;

%   g(7) - velocity gimble (IMU)
%--------------------------------------------------------------------------
g(7)   = x.g.v;

%   g(8:9) - (IMU)
%--------------------------------------------------------------------------
g(8:9) = x.d.x;

return


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
if nargin < 2

    subplot(3,1,1)
    for t = 1:T

        % plot radar
        %------------------------------------------------------------------
        hold off, plot(pP.r.x(1),pP.r.x(2),'^k','MarkerSize',16), hold on

        % plot missile, target and drone (actual)
        %------------------------------------------------------------------
        x = spm_unvec(px(:,t),DEM.G(1).x);
        v = spm_unvec(pv(:,t),DEM.G(2).v);
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        plot(x.m.x(1),x.m.x(2),'.k','MarkerSize',32)
        plot(x.d.x(1),x.d.x(2),'*k','MarkerSize',16)
        plot(x.d.x(1),x.d.x(2),'ok','MarkerSize',8)
        plot(v(1),v(2),'+k','MarkerSize',24)

        x = spm_unvec(qx(:,t),DEM.M(1).x);
        v = spm_unvec(qv(:,t),DEM.M(2).v);
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        % plot missile, target and drone (inferred)
        %------------------------------------------------------------------
        c     = diag(DEM.qU.C{t}) + exp(-4);
        c     = sqrt(c)*3;
        [i,j] = ellipsoid(v(1),v(2),1,c(1),c(2),0,32);

        fill(i(16,:)',j(16,:)',[1 1 1],'LineStyle',':','FaceColor','none','EdgeColor',[1 1 1]/8);
        plot(x.m.x(1),x.m.x(2),'.r','MarkerSize',16)
        plot(x.d.x(1),x.d.x(2),'*r','MarkerSize',8)
        plot(x.d.x(1),x.d.x(2),'or','MarkerSize',4)
        plot(v(1),v(2),'+r','MarkerSize',12)

        % plot gimble
        %------------------------------------------------------------------
        u(1) = x.d.x(1);
        u(2) = x.d.x(1) + sin(x.g.x);
        v(1) = x.d.x(2);
        v(2) = x.d.x(2) + cos(x.g.x);
        plot(u,v,'m')

        % plot bearing
        %------------------------------------------------------------------
        v   = -8;
        plot(g(5),v,'+m','MarkerSize',32)
        plot(g(5),v,'.m','MarkerSize',32)
        plot([-pi,pi],[v,v],':m')

        axis([-1 1 -1 1]*10), axis square
        xlabel('latitude (km)'), ylabel('longitude (km)'), title('Interception')
        drawnow

        % save frame
        %------------------------------------------------------------------
        M(t) = getframe(gca);

    end

    % save movie
    %----------------------------------------------------------------------
    set(gca,'Userdata',[])
    set(gca,'Userdata',{M,8})
    set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

else

    subplot(3,1,2), hold off
    m     = zeros(2,T);
    d     = zeros(2,T);
    for t = 1:T

        % plot missile, target and drone (actual)
        %------------------------------------------------------------------
        x = spm_unvec(px(:,t),DEM.G(1).x);
        v = spm_unvec(pv(:,t),DEM.G(2).v);
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        % distance between missile and drone
        %------------------------------------------------------------------
        m(:,t) = x.m.x;
        d(:,t) = x.d.x;

        hold off, plot(v(1),v(2),'+k','MarkerSize',24), hold on
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
        plot(v(1),v(2),'+r','MarkerSize',12)

        % plot gimble
        %------------------------------------------------------------------
        u(1) = x.d.x(1);
        u(2) = x.d.x(1) + sin(x.g.x)/100;
        v(1) = x.d.x(2);
        v(2) = x.d.x(2) + cos(x.g.x)/100;
        plot(u,v,'m')

        v = DEM.G(2).v;
        axis([-1 1 -1 1]/10 + [v(1) v(1) v(2) v(2)]), axis square
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
    subplot(3,2,5), hold off
    plot(m(1,:),m(2,:)), hold on
    plot(d(1,:),d(2,:)), hold on
    axis([-1 1 -1 1]/100 + [v(1) v(1) v(2) v(2)]), axis square
    xlabel('latitude (km)'), ylabel('longitude (km)'), title('Trajectories (m)')

    % plot distances (after interpolation)
    %----------------------------------------------------------------------
    im(1,:) = interp(m(1,:),1000);
    im(2,:) = interp(m(2,:),1000);
    id(1,:) = interp(d(1,:),1000);
    id(2,:) = interp(d(2,:),1000);

    dx = im - id;
    r  = sqrt(sum(dx.^2));

    subplot(3,2,6), hold off
    semilogy(r*1000,'.'), hold on
    semilogy(r*0 + 1e-0,'-.r')
    semilogy(r*0 + 1e-1,'-.r')
    str = sprintf('Distance (m) (min = %i cm)',fix(min(r)*100000));
    xlabel('time (ms)'), ylabel('meters'), title(str)

end

