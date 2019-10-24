function MDP = MDP_DEM_Oculomotion_Pharma_demo
% Demo of mixed models for oculmotor behaviour, with pharmacological
% manipulations
%
%__________________________________________________________________________
%
% This demo ilustrates the use of mixed (continuous and discrete)
% generative models in simulating oculomotion. An MDP model is used to
% select locations in visual space, and a continuous model is used to
% implement these decisions. See also DEM_demo_MDP_DEM.m,
% MDP_DEM_Oculomotion_demo.m

ACh  = exp(16); % exp(16)
DA   = 1;       % 1
NA   = 5;       % 5
GABA = 1;       % 1  (slow saccades at 1.05)

% Locations
%--------------------------------------------------------------------------
L{1} = [0; -0.2*pi];
L{2} = [0;       0];
L{3} = [0;  0.2*pi];
L{4} = [0.2*pi;  0];
L{5} = [-0.2*pi; 0];

% MDP outcomes to DEM causes
%--------------------------------------------------------------------------
N  = 24;
nl = length(L);
nh = 8;

for i = 1:nh
    for j = 1:nl
        c           = [L{j}; sparse(i,1,1,nh,1)];
        u           = [L{j}; sparse(i,1,1,nh,1)];
        demi.U{i,j} = u*ones(1,N);
        demi.C{i,j} = c*ones(1,N);
    end
end

o     = [6 1];
O{1}  = spm_softmax(sparse(6,1,1,nh,1));
O{2}  = spm_softmax(sparse(2,1,5,nl,1));

% generative model (Continuous)
%==========================================================================
DEM      = DEM_MDP_oculomotion_GM(GABA);
DEM      = spm_MDP_DEM(DEM,demi,O,o);

% generative model (Discrete)
%==========================================================================
mdp      = MDP_DEM_oculomotion_GM(ACh,DA,NA);
mdp.demi = demi;
mdp.DEM  = DEM;

% invert or solve
%--------------------------------------------------------------------------
MDP  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(MDP);

% illustrate belief updating - discrete
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1'); clf
spm_MDP_VB_LFP(MDP,[],2);
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_trial(MDP);

% illustrate belief updating - continuous (last eye movement)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
spm_DEM_qU(MDP.dem(end).qU)

% show oculomotor responses
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4');

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

spm_figure('GetWin','Figure 5'); clf
MDP_plot_delay_task(MDP)

spm_figure('GetWin','Figure 6'); clf
MDP_plot_saccade_characteristics(MDP)

return

function MDP = MDP_DEM_oculomotion_GM(ACh,DA,NA)

T = 4;

% Prior over initial state
%--------------------------------------------------------------------------
D{1} = [0 1 0 0 0]';
D{2} = [1 0 1 1 1]'; % Target                                           
D{3} = [1 0 0]';     % Cue present, fixate, saccade to target           

% Likelihood (Acetylcholine)
%--------------------------------------------------------------------------

for f1 = 1:length(D{1})                                                 
    for f2 = 1:length(D{2})
        for f3 = 1:length(D{3})
            if     f3 == 1
                if f1 == 2
                    A{1}(f2,f1,f2,f3) = 1;  % Cue shown
                else
                    A{1}(7,f1,f2,f3) = 1;   % If not fixating, incorrect
                end
            elseif f3 == 2
                if f1 == 2
                    A{1}(6,f1,f2,f3) = 1;   % No cues shown
                else
                    A{1}(7,f1,f2,f3) = 1;   % If not fixating, incorrect
                end
            else
                if f1 == f2                  
                    A{1}(8,f1,f2,f3) = 1;   % Correct
                else
                    A{1}(7,f1,f2,f3) = 1;   % Incorrect
                end
            end
            A{2}(f1,f1,f2,f3) = 1;
        end
    end
end                                                                     

a{1} = (A{1} + exp(-16)).^ACh;
a{2} = (A{2} + exp(-16)).^ACh;

% Transitions (Noradrenaline)
%--------------------------------------------------------------------------
B{1} = zeros(5);

for k = 1:5
    B{1}(:,:,k) = 0;
    B{1}(k,:,k) = 1;
end
b{1} = B{1};
B{2} = eye(5);
b{2} = (B{2} + exp(-16)).^NA;
B{3} = [0 0 0; 1 0 0; 0 1 1];
b{3} = (B{3} + exp(-16)).^NA;

% Preferences
%--------------------------------------------------------------------------
C{1} = [0 0 0 0 0 0 -6 6]';
C{2} = zeros(5,1);

MDP.A = A;
MDP.a = a;
MDP.B = B;
MDP.b = b;
MDP.C = C;
MDP.D = D;
MDP.V(:,:,1) = [2 2 2 2 2 1 3 4 5;1 2 3 4 5 1 3 4 5;1 2 3 4 5 1 3 4 5];
MDP.V(:,:,2) = ones(size(MDP.V(:,:,1)));
MDP.V(:,:,3) = ones(size(MDP.V(:,:,1)));
MDP.T = T;
MDP.beta = 1/DA;        % Dopamine
MDP.alpha = 0.5;

MDP.tau = 32;
MDP     = spm_MDP_check(MDP);

function DEM = DEM_MDP_oculomotion_GM(GABA)
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
M(1).V  = exp(GABA*[4*ones(1,4) 4*ones(1,4) 16*ones(1,2) 16*ones(1,2) 16*ones(1,5)]);   % error precision (g)
M(1).W  = exp(8);                             % error precision (f)

% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0;0;0;0;0;0;0;0;0;0];                  % priors
M(2).V  = exp(15);


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
G(2).v  = [0;0;0;0;0;0;0;0;0;0];                  % exogenous forces
G(2).a  = spm_vec(a);                         % action forces
G(2).V  = exp(16);


DEM.G = G;
DEM.M = M;

function f = Mf1(x,v,~)
f = [v(1)-x(1); v(2) - x(2); - x(3); - x(4)];

function g = Mg1(x,v,~)
g.l = x;
g.r = x;
g.v = [v(1:2); v(1:2); v(3); v(5); v(6); v(7); v(8)];

function f = Gf1(x,~,a,P)
ke = P.ke;
kv = P.kv;
j  = P.j;
a = spm_unvec(a,P.a);

F.l = a.l - ke*x.l(1:2) - kv*x.l(3:4);
F.r = a.r - ke*x.r(1:2) - kv*x.r(3:4);
f.r = [x.r(3) ; x.r(4); (1/j)*F.r(1); (1/j)*F.r(2)];
f.l = [x.l(3) ; x.l(4); (1/j)*F.l(1); (1/j)*F.l(2)];

function g = Gg1(x,v,~,~)
g.l = x.l;
g.r = x.r;
g.v = [x.l(1:2) ; x.r(1:2); v(3); v(5); v(6); v(7); v(8)];

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

function MDP_plot_delay_task(MDP)
for i = 1:length(MDP.dem)
    dem = MDP.dem(i);
    theta(1,:) = dem.pU.x{1}(1,:);
    theta(2,:) = dem.pU.x{1}(2,:);
    plot(0,0,'+','MarkerSize',50), hold on
    axis([-1 1 -1 1]), axis square, axis off
    for j = 1:size(theta,2)
        if i == 1
            plot(-0.7265,0,'.','MarkerSize',50,'Color',[1 1-dem.Y(13,j) 1-dem.Y(13,j)])
            plot( 0.7265,0,'.','MarkerSize',50,'Color',[1 1-dem.Y(14,j) 1-dem.Y(14,j)])
            plot(0,-0.7265,'.','MarkerSize',50,'Color',[1 1-dem.Y(15,j) 1-dem.Y(15,j)])
            plot( 0,0.7265,'.','MarkerSize',50,'Color',[1 1-dem.Y(16,j) 1-dem.Y(16,j)])
            axis([-1 1 -1 1]), axis square, axis off
        elseif i == 2 && j == 1
            plot(-0.7265,0,'.w','MarkerSize',50)
            plot( 0.7265,0,'.w','MarkerSize',50)
            plot(0,-0.7265,'.w','MarkerSize',50)
            plot( 0,0.7265,'.w','MarkerSize',50)
            axis([-1 1 -1 1]), axis square, axis off
        end
        
        plot(0,0,'+','MarkerSize',50,'Color',[(MDP.s(3,i) == 2) 0 1-(MDP.s(3,i) == 2)])

        plot(tan(theta(2,j)),-tan(theta(1,j)),'.k','MarkerSize',20)
        pause(0.05)
        drawnow
    end
end

function MDP_plot_saccade_characteristics(MDP)
d = [];
for i = 1:length(MDP.dem)
    dem = MDP.dem(i);
    d(end+1:end+size(dem.pU.x{1},2)) = (tan(dem.pU.x{1}(1,:)).^2 + tan(dem.pU.x{1}(2,:)).^2).^0.5;
end
subplot(2,1,1)
plot(d*20/0.7), hold on
title('Distance')
subplot(2,1,2)
plot(length(dem)*gradient(d*20/0.7)*0.25), hold on
title('Speed')
