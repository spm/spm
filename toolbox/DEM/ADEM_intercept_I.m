function ADEM_intercept_I
% Radar tracking and interception in generalised coordinates of motion
%__________________________________________________________________________
% This is a simplified version of the interceptor demo, in which the gimbal
% has a fixed orientation; thereby, reducing the latent states to the
% position of the missile and drone. This enables action to be specified
% directly as controlling velocity (rate of change of position) that is
% projected into generalised coordinates of motion (here to 4th order).
%
% In addition, this demo equips the drone with weaving behaviour modelled
% in terms of a discrete sine set; i.e., parameterised by coefficients that
% control variations in direction of travel as (sine) functions of distance
% to target (i.e., the rendezvous point). It is assumed that these
% coefficients are known or have been previously estimated based upon the
% missiles early behaviour. The remaining latent cause — namely the target
% — is inferred in the usual way; allowing, in principle, en route changes
% to the target location.
%
% The sensor signals are assumed to be sufficiently preprocessed to remove
% random fluctuations that would otherwise confound generalised
% observations. Uncertainty about these signals is therefore modelled in
% terms of the precisions supplied to the generative model. These
% precisions have to be specified carefully to ensure the sensory
% attenuation of observations that can be changed via control. In other
% words, the precision of IMU signals has to be smaller than the precision
% used by the actuator for resolving prediction errors through action.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 2008-2022 Wellcome Centre for Human Neuroimaging

rng(0)

% hidden states, causes and parameters
%==========================================================================
% x    - hidden states:
%  x.m.x    - position (missile) (Km)
%  x.d.x    - position (drone)   (Km)
%
% v    - hidden causes:
%   v  - position (target)
%
% g    - sensor modalities:
%   g(1) - bearing of missile (radar)   (radians)
%   g(2) - log range of missile (radar) (log(Km))
%   g(3) - bearing of missile (IR)      (radians)
%   g(4) - position x (IMU)             (Km)
%   g(5) - position y (IMU)             (Km)
%
% Parameters:
%   P.r.x - position (radar)
%   P.m.s - speed  (missile)
%   P.d.s - speed    (drone)
%--------------------------------------------------------------------------
P.r.x = [4;-2];                               % position  (radar)
P.m.s = 135/60/60;                            % max speed (missile)
P.d.s = 85/60/60;                             % max speed (drone)
P.w   = [1/2;0];                              % evasion weights

% target
%--------------------------------------------------------------------------
v     = -2;                                   % x target location
u     = 0;                                    % prior expecation

% Initial hidden states
%--------------------------------------------------------------------------
x.m.x = [4; 8];                               % position (missile)
x.d.x = [0;-4];                               % position (drone)

% precisions of observation noise
%--------------------------------------------------------------------------
V(1)   = exp(12);                             % bearing (radar)
V(2)   = exp(12);                             % range  (radar)
V(3)   = 0;                                   % bearing (IR)
V(4:5) = exp(2);                              % position (IMU)

% precisions (for action): restricted to IMU signals
%--------------------------------------------------------------------------
U(1)   = 0;                                   % bearing of missile (radar)
U(2)   = 0;                                   % range   of missile (radar)
U(3)   = 0;                                   % bearing of missile (IR gimble)
U(4:5) = exp(8);                              % position (IMU)

% precisions of state noise (assumed by model)
%--------------------------------------------------------------------------
W.m.x = exp([16; 16]);                          % position (missile)
W.d.x = exp([16; 16]);                        % position (drone)

% generative model
%==========================================================================

% order of generalised motion (and smoothness)
%--------------------------------------------------------------------------
dt        = 1;
M(1).E.s  = 1/2;                              % smoothness
M(1).E.n  = 4;                                % order of
M(1).E.d  = 0;                                % generalised motion
M(1).E.dt = dt;                               % time step

% level 1: dynamics of hidden states
%--------------------------------------------------------------------------
M(1).f  = @fx_M;                              % fx (model)
M(1).g  = @gx_M;                              % gx (model)
M(1).x  = x;                                  % hidden states
M(1).V  = spm_vec(V);                         % error precision (signal)
M(1).W  = spm_vec(W);                         % error precision (dynamics)
M(1).pE = P;                                  % parameters

% level 2: hidden cause (rendezvous point - x, y = 0)
%--------------------------------------------------------------------------
M(2).v  = u;                                  % hidden causes
M(2).V  = 0;                                  % hidden causes

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
G(2).v  = v;                                  % hidden cause
G(2).a  = [0; 0];                             % action
G(2).V  = exp(32);


% Assemble DEM structure
%==========================================================================
N       = 240/dt;                             % length of time series (s)
DEM.G   = G;                                  % generative process
DEM.M   = M;                                  % generative model
DEM.C   = repmat(v,1,N);                      % latent causes
DEM.U   = repmat(u,1,N);                      % priors over latent causes

% Solve or integrate
%--------------------------------------------------------------------------
DEM     = spm_ADEM_UNITY(DEM);

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

% update hidden and latent states
%--------------------------------------------------------------------------
DEM.M(1).x = spm_unvec(DEM.qU.x{1}(:,end),DEM.M(1).x);
DEM.M(1).v = spm_unvec(DEM.qU.v{1}(:,end),DEM.M(1).v);
DEM.M(2).v = spm_unvec(DEM.qU.v{2}(:,end),DEM.M(2).v);
DEM.G(1).x = spm_unvec(DEM.pU.x{1}(:,end),DEM.G(1).x);

% increase precision of IR signals
%--------------------------------------------------------------------------
DEM.M(1).V(3,3) = exp(2);                     % precision of IR

% solve
%--------------------------------------------------------------------------
N     = 64/dt;                                % length of time series (s)
DEM.C = repmat(v,1,N);                        % latent causes
DEM.U = repmat(u,1,N);                        % priors over latent causes
DEM   = spm_ADEM_UNITY(DEM);                  % solve

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
%  x.d.x    - position (drone)
%--------------------------------------------------------------------------
f     = x;
x0    = [v(1); 0];                            % target [x,y];
w     = P.w;                                  % DST coeficients;

%  f.m.x    - changes in position (missile)
%--------------------------------------------------------------------------
d     = x0 - x.m.x;
r     = sqrt(d'*d);                           % range
d     = d/r;                                  % direction

% add nonlinear weaving (discrete sine set)
%--------------------------------------------------------------------------
o     = [-d(2); d(1)];
for i = 1:numel(w)
    d = d + w(i)*sin(pi*r*i)*o;
end
r     = sqrt(d'*d);
d     = d/r;
f.m.x = P.m.s*d;

%  f.d.x    - changes in position (drone)
%--------------------------------------------------------------------------
f.d.x = a(1:2);

return


function f = fx_M(x,v,P)
% FORMAT f = fx_M(x,v,P)
% equations of motion for model
%--------------------------------------------------------------------------
%  x.m.x    - position (missile)
%  x.d.x    - position (drone)
%--------------------------------------------------------------------------
f     = x;
x0    = [v(1); 0];                            % target [x,y];
w     = P.w;                                  % DST coeficients;

%  f.m.x    - changes in position (missile)
%--------------------------------------------------------------------------
d     = x0 - x.m.x;
r     = sqrt(d'*d);                           % range
d     = d/r;                                  % direction

% add nonlinear weaving (discrete sine set)
%--------------------------------------------------------------------------
o     = [-d(2); d(1)];
for i = 1:numel(w)
    d = d + w(i)*sin(pi*r*i)*o;
end
r     = sqrt(d'*d);
d     = d/r;
f.m.x = P.m.s*d;

%  f.d.x    - changes in position (drone)
%--------------------------------------------------------------------------
d     = x0 - x.d.x;
r     = sqrt(d'*d);
d     = d/r;
f.d.x = erf(r*8)*P.d.s*d;

return


function g = gx_G(x,v,~,P)
% FORMAT g = gx_G(x,v,a,P)
% observation function (process)
%--------------------------------------------------------------------------
g = gx_M(x,v,P);

return


function g = gx_M(x,~,P)
% FORMAT g = gx_M(x,v,P)
% observation function (model)
%--------------------------------------------------------------------------
%   g(1) - bearing of missile (radar)   (radians)
%   g(2) - log range of missile (radar) (log(Km))
%   g(3) - bearing of missile (IR)      (radians)
%   g(4) - position x (IMU)             (Km)
%   g(5) - position y (IMU)             (Km)
%--------------------------------------------------------------------------
g      = zeros(5,1);

%   g(1) - bearing of missile (radar)
%--------------------------------------------------------------------------
dx     = x.m.x - P.r.x;
g(1)   = atan(dx(1)/dx(2));

%   g(2) - range of missile (radar)
%--------------------------------------------------------------------------
r      = sqrt(dx'*dx);
g(2)   = log(r);

%   g(3) - bearing of missile (IR gimble)
%--------------------------------------------------------------------------
dx     = x.m.x - x.d.x;
g(3)   = atan(dx(1)/dx(2));

%   g(4:5) - position (IMU) 
%--------------------------------------------------------------------------
g(4:5) = x.d.x;


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

    subplot(2,2,1)
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
        plot(v(1),0,'+k','MarkerSize',24)

        x = spm_unvec(qx(:,t),DEM.M(1).x);
        v = spm_unvec(qv(:,t),DEM.M(2).v);
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        % plot missile, target and drone (inferred)
        %------------------------------------------------------------------
        c     = DEM.qU.C{t}(1);
        c     = sqrt(c)*3;
        [i,j] = ellipsoid(v(1),0,1,c(1),1/4,0,32);

        fill(i(16,:)',j(16,:)',[1 1 1],'LineStyle',':','FaceColor','none','EdgeColor',[1 1 1]/8);
        plot(x.m.x(1),x.m.x(2),'.r','MarkerSize',16)
        plot(x.d.x(1),x.d.x(2),'*r','MarkerSize',8)
        plot(x.d.x(1),x.d.x(2),'or','MarkerSize',4)
        plot(v(1),0,'+r','MarkerSize',12)

        % plot gimble
        %------------------------------------------------------------------
        u(1) = x.d.x(1);
        u(2) = x.d.x(1) + sin(0);
        v(1) = x.d.x(2);
        v(2) = x.d.x(2) + cos(0);
        plot(u(1:2),v(1:2),'m')

        % plot bearing
        %------------------------------------------------------------------
        v   = -8;
        plot(g(3)*8,v,'+m','MarkerSize',32)
        plot(g(3)*8,v,'.m','MarkerSize',32)
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
        v = spm_unvec(pv(:,t),DEM.G(2).v);
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
        plot(v(1),0,'+r','MarkerSize',12)

        % plot gimble
        %------------------------------------------------------------------
        u(1) = x.d.x(1);
        u(2) = x.d.x(1) + sin(0)/100;
        v(1) = x.d.x(2);
        v(2) = x.d.x(2) + cos(1)/100;
        plot(u,v,'m')

        v = DEM.G(2).v;
        axis([-1 1 -1 1] + [v(1) v(1) 0 0]), axis square
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
    plot(m(1,:),m(2,:),'--k'), hold on
    plot(d(1,:),d(2,:),'--r'), hold on
    axis([-1 1 -1 1] + [v(1) v(1) 0 0]), axis square
    xlabel('latitude (km)'), ylabel('longitude (km)'), title('Trajectories (m)')

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

