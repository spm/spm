function ADEM_intercept_III
% Fixed wing interception in generalised coordinates of motion
%__________________________________________________________________________
% This version is for a fixed-wing interceptor. In brief, it uses a
% combination of first order kinetics to approach the bomb in a moving
% frame of reference. This moving frame of reference is modelled using the
% second order dynamics based upon a hidden cause that pulls both the
% interceptor and bomb to a common point, which is moving from towards the
% target. This hidden causes equipped with smooth random fluctuations to
% model countermeasures.
%__________________________________________________________________________

% Copyright (C) Karl Friston

rng(1)

% hidden states, causes and parameters
%==========================================================================
% x    - hidden states:
%  x.m.x    - position (missile) (Km)
%  x.m.v    - velocity (missile) (Km/s)
%  x.d.x    - position (drone)   (Km)
%  x.d.v    - velocity (drone)   (Km/s)
%
% v    - hidden causes: x,y, position of an attracting point
%
% g    - sensor modalities (for one radar):
%   g(:) - bearing of missile   (radar) (radians)
%   g(:) - elevation of missile (radar) (Km)
%   g(:) - log range of missile (radar) (log(Km))
%   g(:) - position x (IMU)             (Km)
%   g(:) - position y (IMU)             (Km)
%   g(:) - position z (IMU)             (Km)
%   g(:) - yaw        (IMU)             (rad)
%   g(:) - pitch      (IMU)             (rad)
%   g(:) - range x (IR)                 (Km)
%   g(:) - range y (IR)                 (Km)
%   g(:) - range z (IR)                 (Km)
%   g(:) - range of missile (IR)        (Km)
%--------------------------------------------------------------------------

% Generative model
%==========================================================================
% First, specify the parameters of the missile and drone
%--------------------------------------------------------------------------

% parameters
%--------------------------------------------------------------------------
P.r.x = [4, -4;-2, -1; 0, 0];                 % position (radars)
P.d.s = 400/60/60;                            % max speed (drone) per sec
P.m.s = 400/60/60;                            % max speed (missile) per sec
P.m.m = 1/16;                                 % mass (missile)
P.d.m = 1/16;                                 % mass (missile)
P.m.r = 2;                                    % pull range (missile) (km)
P.d.r = 2;                                    % pull range (drone) (km)

% hidden states (initial)
%--------------------------------------------------------------------------
x.m.x = [4;12; 2];                            % position (missile)
x.m.v = [0; 0; 0];                            % velocity (missile)
x.d.x = [0; 0; 0];                            % position (drone)
x.d.v = [0; 0; 0];                            % velocity (missile)

% precisions of state noise (assumed by model)
%--------------------------------------------------------------------------
W.m.x = exp([16; 16; 16]);                    % position (missile)
W.m.v = exp([16; 16; 16]);                    % velocity (missile)
W.d.x = exp([16; 16; 16]);                    % position (drone)
W.d.v = exp([16; 16; 16]);                    % velocity (drone)

% precisions of observation noise
%--------------------------------------------------------------------------
n     = size(P.r.x,2);                        % number of radar
V.rad = zeros(3,n) + exp(12);                 % range & bearing (radar)
V.imu = zeros(3,1) + exp(-2);                 % drone (IMU)
V.yaw = zeros(2,1);                           % drone (IMU)
V.b   = zeros(3,1) + exp(4);                  % position (IR)
V.r   = 0;                                    % range (IR)

% precisions (for action): restricted to IMU signals
%--------------------------------------------------------------------------
U     = spm_zeros(V);                         % action precision
U.imu = U.imu + exp(4);                       % position (IMU)

% generative model for ADEM
%==========================================================================

% order of generalised motion (and smoothness)
%--------------------------------------------------------------------------
dt     = 1/10;                                % time step (seconds)
E.s    = 1/2;                                 % smoothness (seconds)
E.n    = 3;                                   % order (states)
E.d    = 0;                                   % order (causes)
E.dt   = dt;                                  % time step

% level 1: dynamics of hidden states
%--------------------------------------------------------------------------
M(1).f  = @fx_M;                              % fx (model)
M(1).g  = @gx_M;                              % gx (model)
M(1).x  = x;                                  % hidden states
M(1).V  = spm_vec(V);                         % error precision (signal)
M(1).W  = spm_vec(W);                         % error precision (dynamics)
M(1).pE = P;                                  % parameters
M(1).E  = E;

% level 2: hidden cause
%--------------------------------------------------------------------------
M(2).v  = [0; 0; 0];                          % hidden causes
M(2).V  = 0;                                  % hidden causes

% generative process
%==========================================================================

% Trajectory of missile
%--------------------------------------------------------------------------
C0      = x.m.x - [0; P.m.r; 0];              % start (ahead of missile)
CT      = [6; 4; 1];                          % stop (target)
d       = CT - C0;                            % distance
N       = sqrt(d'*d)/P.m.s/dt;                % length of time series (s)
for i = 1:numel(C0)
    C(i,:) = linspace(C0(i),CT(i),N);
end

% add smooth random fluctuations
%--------------------------------------------------------------------------
R       = randn(size(C));
R       = spm_conv(R,0,4/dt)/dt;
C       = C + R/2;

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
G(2).v  = C(:,1);                             % hidden causes
G(2).a  = [0; 0; 0];                          % action
G(2).V  = exp(32);

% Initial interception
%==========================================================================
% Having specified the generative model, and process, we can now simulate an
% initial interception that pauses when the interceptor is within 1 km of
% the bomb
%--------------------------------------------------------------------------

% Assemble DEM structure
%--------------------------------------------------------------------------
DEM.G   = G;                                  % generative process
DEM.M   = M;                                  % generative model
DEM.C   = C;                                  % latent causes
DEM.U   = spm_zeros(C);                       % priors over latent causes

% Solve or integrate
%--------------------------------------------------------------------------
DEM.break = 'qu.y{1}(end)  < 1';              % break if within 1 km
DEM       = spm_ADEM_UNITY(DEM);

% show results
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
% We can now resume the simulation with IR information now available,
% because the bone is sufficiently close to the interceptor. We will
% terminate simulation when the interceptor is within 1 m or the target is
% reached
%--------------------------------------------------------------------------
DEM.break  = 'qu.y{1}(end) < 1e-3';           % break if within 1 m

% update hidden states and causes
%--------------------------------------------------------------------------
DEM.M(1).x = spm_unvec(DEM.qU.x{1}(:,end),DEM.M(1).x);
DEM.M(1).v = spm_unvec(DEM.qU.v{1}(:,end),DEM.M(1).v);
DEM.M(2).v = spm_unvec(DEM.qU.v{2}(:,end),DEM.M(2).v);
DEM.G(1).x = spm_unvec(DEM.pU.x{1}(:,end),DEM.G(1).x);
DEM.G(2).v = spm_unvec(DEM.pU.v{2}(:,end),DEM.G(2).v);

% solve
%--------------------------------------------------------------------------
n        = size(DEM.Y,2);
C(:,1:n) = [];                               % remove past causes
DEM.C    = C;                                % latent causes
DEM.U    = spm_zeros(C);                     % priors over latent causes
FEM      = spm_ADEM_UNITY(DEM);              % solve

% create movie of radar tracking and final approach
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_graphics_intercept(FEM,1)

return


% model functions (subroutines)
%__________________________________________________________________________

function f = fx_G(x,v,a,P)
% FORMAT f = fx_G(x,v,a,P)
% equations of motion for generative process
%--------------------------------------------------------------------------
%  x.m.x    - position (missile)
%  x.d.x    - position (drone)
%--------------------------------------------------------------------------
f     = fx_missile(x,v,P);

% f.d.x - changes in position (interceptor)
%--------------------------------------------------------------------------
f.d.x = erf(a);

return


function f = fx_missile(x,v,P)
% FORMAT f = fx_M(x,v,P)
% equations of motion for missile
% These are second-order dynamics in which the missile is effectively
% attracted to a fixed point (hidden cause) with a spring whose stiffness
% is calibrated to ensure a steady-state solution in which the missile lags
% the (moving) fixed point by P.m.r
%--------------------------------------------------------------------------
%  x.m.x - position (missile)
%  x.d.x - position (drone)
%--------------------------------------------------------------------------
f     = x;

%  f.m.x - changes in position (missile)
%--------------------------------------------------------------------------
f.m.x = x.m.v;

% f.m.x - changes in velocity (missile)
%--------------------------------------------------------------------------
d     = v - x.m.x;
f.m.v = (P.m.s*d/P.m.r - x.m.v)/P.m.m;


function f = fx_M(x,v,P)
% FORMAT f = fx_M(x,v,P)
% equations of motion for model
% This is the same as the missile with an additional first-order change in
% velocity designed to approach the missile at 50 km/h, slowing down when
% within 100 m.
%--------------------------------------------------------------------------
%  x.m.x - position (missile)
%  x.d.x - position (drone)
%--------------------------------------------------------------------------
f     = fx_missile(x,v,P);

% f.d.x - changes in position (interceptor)
%--------------------------------------------------------------------------
d     = x.m.x - x.d.x;
r     = sqrt(d'*d);
d     = d/r;
f.d.x = x.d.v + 50e-3*erf(r/100e-3)*d;

% f.d.v - changes in velocity (interceptor)
%--------------------------------------------------------------------------
d     = v - x.d.x;
f.d.v = (P.m.s*d/P.d.r - x.d.v)/P.d.m;

return


function g = gx_G(x,v,a,P)
% FORMAT g = gx_G(x,v,a,P)
% observation function (process)
%--------------------------------------------------------------------------
x.d.v = erf(a);
g     = gx_M(x,v,P);

return

function g = gx_M(x,v,P)
% FORMAT g = gx_M(x,v,P)
% observation function (model)
% When the missile is out of range, the IR just reports a relative position
% of zero, and a fixed range of 8 km. The partial derivatives of this
% output with respect to causes and states are therefore zero and they do
% not contribute to control as inference
%--------------------------------------------------------------------------
%   g(1) - bearing of missile (radar)   (radians)
%   g(2) - log range of missile (radar) (log(Km))
%   g(3) - position x (IMU)             (Km)
%   g(4) - position y (IMU)             (Km)
%   ...
%--------------------------------------------------------------------------
g   = gx_radar(x,v,P);

%   g - position (IMU)
%--------------------------------------------------------------------------
g(end + 1) = x.d.x(1);
g(end + 1) = x.d.x(2);
g(end + 1) = x.d.x(3);

%   g - yaw and pitch (IMU)
%--------------------------------------------------------------------------
dx    = x.d.v;
g(end + 1) = atan2(dx(1),dx(2));
g(end + 1) = asin(dx(3)/sqrt(sum(dx.^2) + eps));


%   g - bearing and range of missile (IR gimble)
%--------------------------------------------------------------------------
d  = x.m.x - x.d.x;
r  = sqrt(d'*d);
if r < 1
    g(end + 1) = d(1);
    g(end + 1) = d(2);
    g(end + 1) = d(3);
    g(end + 1) = r;
else
    g(end + 1) = 0;
    g(end + 1) = 0;
    g(end + 1) = 0;
    g(end + 1) = 8;
end

return

function g = gx_radar(x,v,P)
% FORMAT g = gx_O(x,v,P)
% observation function (radar)
%--------------------------------------------------------------------------
%   g(1,i) - bearing of missile (radar)   (radians)
%   g(2,i) - log range of missile (radar) (log(Km))
%--------------------------------------------------------------------------
n     = size(P.r.x,2);
g     = zeros(2,n);
for i = 1:n

    %   g(i,1) - bearing of missile (radar)
    %--------------------------------------------------------------------------
    d      = x.m.x - P.r.x(:,i);
    g(1,i) = atan2(d(1),d(2));

    %   g(i,2) - elevation (radar)
    %--------------------------------------------------------------------------
    g(2,i) = x.m.x(3);


    %   g(i,3) - range of missile (radar)
    %--------------------------------------------------------------------------
    r      = sqrt(d'*d);
    g(3,i) = log(r);

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
iy = (1:2) + nr*3 + 3;
if nargin < 2

    subplot(2,1,1)
    for t = 1:T

        % plot radar
        %------------------------------------------------------------------
        hold off,
        for i = 1:nr
            plot(pP.r.x(1,i),pP.r.x(2,i),'^k','MarkerSize',16), hold on
        end

        % plot target
        %------------------------------------------------------------------
        plot(DEM.C(1,end),DEM.C(2,end),'*g','MarkerSize',32), hold on

        % plot missile, target and drone (actual)
        %------------------------------------------------------------------
        x = spm_unvec(px(:,t),DEM.G(1).x);
        v = spm_unvec(pv(:,t),DEM.G(2).v);
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        % missile and drone positions
        %------------------------------------------------------------------
        m(:,t) = x.m.x;
        d(:,t) = x.d.x;

        plot(x.m.x(1),x.m.x(2),'.k','MarkerSize',32)
        plot(x.d.x(1),x.d.x(2),'ok','MarkerSize',16)
        plot(v(1),v(2),'+k','MarkerSize',24)

        x = spm_unvec(qx(:,t),DEM.M(1).x);
        v = spm_unvec(qv(:,t),DEM.M(2).v);

        % plot missile, target and drone (inferred)
        %------------------------------------------------------------------
        c     = DEM.qU.C{t};
        c     = sqrt(c)*3;
        [i,j] = ellipsoid(v(1),v(2),1,c(1,1),c(2,2),0,32);

        fill(i(16,:)',j(16,:)',[1 1 1],'LineStyle',':','FaceColor','none','EdgeColor',[1 1 1]/8);
        plot(x.m.x(1),x.m.x(2),'.r','MarkerSize',16)
        plot(x.d.x(1),x.d.x(2),'or','MarkerSize',8)
        plot(v(1),v(2),'+r','MarkerSize',12)

        % plot direction
        %------------------------------------------------------------------
        a    = atan2(x.d.v(1),x.d.v(2));
        p(1) = x.d.x(1);
        p(2) = x.d.x(1) + sin(a);
        q(1) = x.d.x(2);
        q(2) = x.d.x(2) + cos(a);
        plot(p,q,'r')

        a    = atan2(x.m.v(1),x.m.v(2));
        p(1) = x.m.x(1);
        p(2) = x.m.x(1) + sin(a);
        q(1) = x.m.x(2);
        q(2) = x.m.x(2) + cos(a);
        plot(p,q,'r')

        % plot bearing
        %------------------------------------------------------------------
        b   = atan2(g(end - 3),g(end - 2));
        v   = -6;
        plot(b,v,'+m','MarkerSize',32)
        plot(b,v,'.m','MarkerSize',32)
        plot([-pi,pi],[v,v],':m')

        axis([-1 1 -1 1]*12), axis square
        xlabel('latitude (km)'), ylabel('longitude (km)'),
        title('Interception')

        % save frame
        %------------------------------------------------------------------
        drawnow
        % M(t) = getframe(gca);

    end

    % save movie
    %----------------------------------------------------------------------
    % set(gca,'Userdata',[])
    % set(gca,'Userdata',{M,8})
    % set(gca,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

    % plot trajectories
    %----------------------------------------------------------------------
    plot(m(1,:),m(2,:),':k'), hold on
    plot(d(1,:),d(2,:),':r'), hold on
    drawnow

else

    m     = zeros(3,T);
    d     = zeros(3,T);
    for t = 1:T

        % plot missile, target and drone (actual)
        %------------------------------------------------------------------
        x = spm_unvec(px(:,t),DEM.G(1).x);
        v = spm_unvec(pv(:,t),DEM.G(2).v);
        g = spm_unvec(gv(:,t),DEM.G(1).v);

        % missile
        %------------------------------------------------------------------
        m(:,t) = x.m.x;
        d(:,t) = x.d.x;

        subplot(2,2,1), hold off
        plot(v(1),v(2),'+k','MarkerSize',24), hold on
        plot(x.m.x(1),x.m.x(2),'.k','MarkerSize',64)
        plot(x.d.x(1),x.d.x(2),'ok','MarkerSize',16)

        % target
        %------------------------------------------------------------------
        plot(DEM.C(1,end),DEM.C(2,end),'*g','MarkerSize',16), hold on

        % plot direction
        %------------------------------------------------------------------
        a    = g(iy(1));
        p(1) = x.d.x(1);
        p(2) = x.d.x(1) + sin(a)/2;
        q(1) = x.d.x(2);
        q(2) = x.d.x(2) + cos(a)/2;
        plot(p,q,'k')

        a    = atan2(x.m.v(1),x.m.v(2));
        p(1) = x.m.x(1);
        p(2) = x.m.x(1) + sin(a)/2;
        q(1) = x.m.x(2);
        q(2) = x.m.x(2) + cos(a)/2;
        plot(p,q,'k')

        % and drone positions
        %------------------------------------------------------------------
        x = spm_unvec(qx(:,t),DEM.M(1).x);
        v = spm_unvec(qv(:,t),DEM.M(2).v);

        % plot missile, target and drone (inferred)
        %------------------------------------------------------------------
        plot(v(1),v(2),'+r','MarkerSize',24)
        plot(x.m.x(1),x.m.x(2),'.r','MarkerSize',32)
        plot(x.d.x(1),x.d.x(2),'or','MarkerSize',16)

        % plot direction
        %------------------------------------------------------------------
        a    = atan2(x.d.v(1),x.d.v(2));
        p(1) = x.d.x(1);
        p(2) = x.d.x(1) + sin(a)/2;
        q(1) = x.d.x(2);
        q(2) = x.d.x(2) + cos(a)/2;
        plot(p,q,'r')

        a    = atan2(x.m.v(1),x.m.v(2));
        p(1) = x.m.x(1);
        p(2) = x.m.x(1) + sin(a)/2;
        q(1) = x.m.x(2);
        q(2) = x.m.x(2) + cos(a)/2;
        plot(p,q,'r')

        axis([-1 1 -1 1]*3 + [v(1) v(1) v(2) v(2)]), axis square, grid on
        xlabel('latitude (km)'), ylabel('longitude (km)'), title('Interception')


        % plot camera views
        %==================================================================

        % process states: relative bearing and elevation in radians
        %------------------------------------------------------------------
        x  = spm_unvec(px(:,t),DEM.G(1).x);
        dx = x.m.x - x.d.x;
        db = atan2(dx(1),dx(2));
        de = asin(dx(3)/sqrt(sum(dx(1:2).^2)) + eps);

        rb = db - g(iy(1));
        if rb > +pi, rb = rb - 2*pi; end
        if rb < -pi, rb = rb + 2*pi; end
        re = de - g(iy(2));
        rr = sqrt(dx'*dx);
        rs = 4*atan(1e-3/rr);
        ri = [-1  1  1 -1];
        rj = [-1 -1  1  1];


        % front view
        %------------------------------------------------------------------
        subplot(4,6,5)
        fill(rb + ri*rs,re + rj*rs,'k')
        axis([-1 1 -1 1]*pi/2), axis square, title('front view')

        % left view
        %------------------------------------------------------------------
        subplot(4,6,4)
        fill(rb + ri*rs,re + rj*rs,'k')
        axis([-2 0 -1 1]*pi/2), axis square, title('left')

        % right view
        %------------------------------------------------------------------
        subplot(4,6,6)
        fill(rb + ri*rs,re + rj*rs,'k')
        axis([0 2 -1 1]*pi/2), axis square, title('right')

        % rear view
        %------------------------------------------------------------------
        subplot(4,6,6 + 5)
        if rb > 0, rb = rb - pi; end
        if rb < 0, rb = rb + pi; end
        fill(rb + ri*rs,re + rj*rs,'k')
        axis([-1 1 -1 1]*pi/2), axis square, title('rear view')

        % save frame
        %------------------------------------------------------------------
        drawnow
        % M(t) = getframe(gcf);

    end

    % save movie
    %----------------------------------------------------------------------
    %     set(gcf,'Userdata',[])
    %     set(gcf,'Userdata',{M,8})
    %     set(gcf,'ButtonDownFcn','spm_DEM_ButtonDownFcn')

    % plot speeds
    %----------------------------------------------------------------------
    subplot(4,3,7), hold off
    dt = DEM.M(1).E.dt;
    t  = (1:T)*dt;
    vm = gradient(m,dt);
    vm = sqrt(sum(vm.^2));
    vd = gradient(d,dt);
    vd = sqrt(sum(vd.^2));
    plot(t/60,vm*3600,'k'), hold on
    plot(t/60,vd*3600,'r'), hold on
    axis square
    xlabel('time (min)'), ylabel('Km/h'), title('Speeds')


    % plot trajectories
    %----------------------------------------------------------------------
    subplot(4,3,8), hold off
    plot(m(1,:),m(2,:),'k'), hold on
    plot(d(1,:),d(2,:),'r'), hold on
    axis square
    xlabel('latitude (km)'), ylabel('longitude (km)'), title('Trajectories')

    % plot distances (after interpolation)
    %----------------------------------------------------------------------
    im(1,:) = interp(m(1,:),1000);
    im(2,:) = interp(m(2,:),1000);
    id(1,:) = interp(d(1,:),1000);
    id(2,:) = interp(d(2,:),1000);

    dx = im - id;
    r  = sqrt(sum(dx.^2));

    subplot(4,3,9), hold off
    semilogy(r*1000,'.'), hold on
    semilogy(r*0 + 1e-0,'-.r')
    semilogy(r*0 + 1e-1,'-.r')
    str = sprintf('Distance (m) (min = %.2f m)',min(r)*1000);
    xlabel('time (ms)'), ylabel('meters'), title(str)
    axis square

end
