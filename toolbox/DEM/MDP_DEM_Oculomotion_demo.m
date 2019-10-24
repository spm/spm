function MDP = MDP_DEM_Oculomotion_demo
% Demo of mixed models for oculmotor behaviour
%
%__________________________________________________________________________
%
% This demo ilustrates the use of mixed (continuous and discrete)
% generative models in simulating oculomotion. An MDP model is used to
% select locations in visual space, and a continuous model is used to
% implement these decisions. See also DEM_demo_MDP_DEM.m.
% For a version of this routine with simulated pharmacological
% interventions (and a delay-period task) please see: 
% MDP_DEM_Oculomotion_Pharma_demo.m

% Locations
%--------------------------------------------------------------------------
L{1} = [0; -0.2*pi];
L{2} = [0;       0];
L{3} = [0;  0.2*pi];

% MDP outcomes to DEM causes
%--------------------------------------------------------------------------
N  = 24;
nl = length(L);
nh = 2;

for i = 1:nh
    for j = 1:nl
        c           = [L{j}; sparse(i,1,1,nh,1)];
        u           = [L{j}; sparse(i,1,1,nh,1)];
        demi.U{i,j} = u*ones(1,N);
        demi.C{i,j} = c*ones(1,N);
    end
end

o     = [1 1];
O{1}  = spm_softmax(sparse(1:2,1,1,nh,1));
O{2}  = spm_softmax(sparse(3,1,3,nl,1));

% generative model (Continuous)
%==========================================================================
DEM      = DEM_MDP_oculomotion_GM;
DEM      = spm_MDP_DEM(DEM,demi,O,o);

% generative model (Discrete)
%==========================================================================
mdp      = MDP_DEM_oculomotion_GM;
mdp.demi = demi;
mdp.DEM  = DEM;

% invert or solve
%--------------------------------------------------------------------------
MDP  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(MDP);

% illustrate belief updating - discrete
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_LFP(MDP);

% illustrate belief updating - continuous (last eye movement)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf;
spm_DEM_qU(MDP.dem(end).qU)

% show oculomotor responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');
for i = 1:length(MDP.dem)
    dem = MDP.dem(i);
    thetal(1,:) = dem.pU.x{1}(1,:);
    thetal(2,:) = dem.pU.x{1}(2,:);
    thetal(3,:) = zeros(1,length(thetal));
    thetar(1,:) = dem.pU.x{1}(5,:);
    thetar(2,:) = dem.pU.x{1}(6,:);
    thetar(3,:) = zeros(1,length(thetar));
    
    subplot(2,1,1)
    MDP_plot_eyes(thetal,thetar);
    MDP_plot_collicular_map(full(dem.qU.z{2}))
end


return

function MDP = MDP_DEM_oculomotion_GM

T = 4;

% Prior over initial state
%--------------------------------------------------------------------------
D{1} = [0,1,0]';

% Likelihood
%--------------------------------------------------------------------------
A{1} = ones(2,3);
A{2} = eye(3);

% Transitions
%--------------------------------------------------------------------------
B{1} = zeros(3);

for k = 1:3
    B{1}(:,:,k) = 0;
    B{1}(k,:,k) = 1;
end

% Preferences
%--------------------------------------------------------------------------
C{1} = zeros(2,T);
C{2} = [0 0 1 0;
    1 0 0 1;
    0 1 0 0];

MDP.A = A;
MDP.B = B;
MDP.C = C;
MDP.D = D;
MDP.T = T;

MDP.tau = 32;
MDP     = spm_MDP_check(MDP);

function DEM = DEM_MDP_oculomotion_GM
% hidden states
%--------------------------------------------------------------------------
x.l = [0;0;0;0];       % left eye oculomotor angle (v;h) and velocity (v;h)
x.r = [0;0;0;0];       % right eye

% and actions
%--------------------------------------------------------------------------
a.l = [0;0];
a.r = [0;0];

% parameters
%--------------------------------------------------------------------------
P.ke = 2;              % elasticty of tendons
P.kv = 1;              % viscocity of orbit
P.a = a;
P.j = 1.5;

M(1).pE = P;

% recognition model
%==========================================================================
M(1).E.s = 1/2;                               % smoothness: default 1/2
M(1).E.n = 4;                                 % order of
M(1).E.d = 2;                                 % generalised motion

% level 1
%--------------------------------------------------------------------------
M(1).f  = @(x,v,P) Mf1(x,v,P);                % plant dynamics
M(1).g  = @(x,v,P) Mg1(x,v,P);                % prediction

M(1).x  = x.l;                                % hidden states
M(1).V  = exp([4*ones(1,4) 4*ones(1,4) 16*ones(1,2) 16*ones(1,2)]);   % error precision (g)
M(1).W  = exp(8);                             % error precision (f)

% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0;0;0;0];                          % priors
M(2).V  = exp(16);


% generative model
%==========================================================================

% first level
%--------------------------------------------------------------------------
G(1).f  = @(x,v,a,P) Gf1(x,v,a,P);
G(1).g  = @(x,v,a,P) Gg1(x,v,a,P);
G(1).x  = x;                                  % hidden states
G(1).V  = exp(16);                            % error precision
G(1).W  = exp(16);                            % error precision
G(1).U  = exp(8);                             % gain
G(1).pE = P;

% second level
%--------------------------------------------------------------------------
G(2).v  = [0;0;0;0];                          % exogenous forces
G(2).a  = spm_vec(a);                         % action forces
G(2).V  = exp(16);


DEM.G = G;
DEM.M = M;

function f = Mf1(x,v,~)
f = [v(1)-x(1); v(2) - x(2); - x(3); - x(4)];

function g = Mg1(x,v,~)
g.l = x;
g.r = x;
g.v = [v(1:2); v(1:2)];

function f = Gf1(x,~,a,P)
ke = P.ke;
kv = P.kv;
j  = P.j;
a = spm_unvec(a,P.a);

F.l = a.l - ke*x.l(1:2) - kv*x.l(3:4);
F.r = a.r - ke*x.r(1:2) - kv*x.r(3:4);
f.r = [x.r(3) ; x.r(4); (1/j)*F.r(1); (1/j)*F.r(2)];
f.l = [x.l(3) ; x.l(4); (1/j)*F.l(1); (1/j)*F.l(2)];

function g = Gg1(x,~,~,~)
g.l = x.l;
g.r = x.r;
g.v = [x.l(1:2) ; x.r(1:2)];

% Plotting routines
%--------------------------------------------------------------------------
function MDP_plot_eyes(thetal,thetar)


[x,y,z] = sphere;
surf(x,y,z,'Facecolor','w','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7)
hold on
surf(x + 3,y,z,'Facecolor','w','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7)

a = hgtransform;
b = hgtransform;

xi = (x(1:5,:)); yi = y(1:5,:); zi = z(1:5,:)-exp(-4);
surf(xi,yi,zi,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',a)
surf(xi,yi,zi,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',b)

xp = x(1:3,:); yp = y(1:3,:); zp = z(1:3,:)-exp(-4);
surf(xp,yp,zp,'Facecolor','k','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',a)
surf(xp,yp,zp,'Facecolor','k','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',b)

light
view(2)
axis equal
axis off

for i = 1:size(thetal,2)
    
    Rxa = makehgtform('xrotate',pi+thetal(1,i));
    Rxb = makehgtform('xrotate',pi+thetar(1,i));
    Rya = makehgtform('yrotate',thetal(2,i));
    Ryb = makehgtform('yrotate',thetar(2,i));
    Rza = makehgtform('zrotate',thetal(3,i));
    Rzb = makehgtform('zrotate',thetar(3,i));
    Rt = makehgtform('translate', [3 0 0]);
    
    
    set(a,'Matrix',Rxa*Rya*Rza);
    set(b,'Matrix',Rt*Rxb*Ryb*Rzb);
    drawnow
    pause(0.05);
    
end
hold off

function MDP_plot_collicular_map(theta)
H = (theta(2,:));
V = (theta(1,:));

x = 0:pi/16:pi/2;
y = -pi/2:pi/4:pi/2;
xf = fliplr(x) - pi/2;
for j = 1:length(theta)
    
    subplot(2,2,4)
    for i = 1:length(y)
        plot(x,y(i)*log(50*x+1),'--k'), hold on
    end
    for i = 1:length(x)
        plot(x(i)*ones(size(y)),y*log(50*x(i)+1),'--k')
    end
    axis square
    axis off
    subplot(2,2,3)
    for i = 1:length(y)
        plot(xf,y(i)*log(50*x+1),'--k'), hold on
    end
    for i = 1:length(x)
        plot(xf(i)*ones(size(y)),y*log(50*x(i)+1),'--k')
    end
    axis square
    axis off
    if H(j)>0
        subplot(2,2,4)
    else
        subplot(2,2,3)
    end
    plot(H(j),V(j)*log(50*H(j)+1),'.r','MarkerSize',50), hold off
    axis square
    axis off
    pause(0.05)
end
