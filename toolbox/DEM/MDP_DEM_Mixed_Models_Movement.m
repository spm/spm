function MDP_DEM_Mixed_Models_Movement
% This demo illustrates a series of computational pathologies as elicited
% during a synthetic neurological examination. This focuses upon an
% examination of the biceps reflex, and a simple coordination task. Each of
% these are simulated for an arm that may move in three dimensions through
% internal and external rotation of the shoulder, flexion and extension of
% the shoulder, and flexion and extension of the elbow. These dynamics play
% out through an active Bayesian filtering scheme, where a series of
% attracting points draw the hand to different locations in 3D space. The
% selection of these attracting points involves a hierarhical Markov
% Decision Process, which identifies these sequences based upon the prior
% belief that (1) the sequence will minimise expected free energy and (2)
% the sequence is consistent with the trajectory predicted by the highest
% (contextual) level of the model.
%__________________________________________________________________________


rng default

% REFLEXES
%==========================================================================
% Demonstration 1 - Tendon reflexes (healthy, pyramidal lesion,
% cerebellar lesion)

% Specify generative model and process
%--------------------------------------------------------------------------

a = 2;                                      % length of upper arm
b = 1.5;                                    % length of forearm
c = [-1;-1;-1];                             % location of shoulder

DEM        = DEM_MDP_Movement_Disorders_GM; % Continuous state generative model
x = [0;-1;1;0;0;0];                         % Initial joint positions


DEM.G(1).g = @(x,v,a,P) Gg1Reflex(x,v,a,P); % Generative process allows for tendon tap
DEM.M(1).x = x;
DEM.G(1).x = x;
DEM.U = [MDP_DEM_MD_transform(x(1:3),a,b,c);0;0;0]*ones(1,64);
DEM.C = zeros(size(DEM.U));
DEM.C(3,10:12) = -0.25;                     % Transient stimulation of type II sensory afferents

% Solve model and plot healthy reflexes
%--------------------------------------------------------------------------
dem = spm_ADEM(DEM);

% Animate
spm_figure('GetWin','Figure 1'); clf
opengl('software')
for i = 1:size(dem.pU.x{1},2)
    MDP_DEM_MD_plot_arm(full(dem.pU.x{1}(1:3,i)),dem.M(1).pE.b1,dem.M(1).pE.b2,[-1;-1;-1]), hold on
    light
    view(10,0)
    zlim([-4 1])
    xlim([-1.2 2.1])
    drawnow
    hold off
    title('Normal')
end

clear dem

% Upper motor lesion
%--------------------------------------------------------------------------
DEM.M(1).V = exp(5);        % Overestimate precision

dem = spm_ADEM(DEM);

% Animate
spm_figure('GetWin','Figure 1'); clf
for i = 1:size(dem.pU.x{1},2)
    MDP_DEM_MD_plot_arm(full(dem.pU.x{1}(1:3,i)),dem.M(1).pE.b1,dem.M(1).pE.b2,[-1;-1;-1]), hold on
    light
    view(10,0)
    zlim([-4 1])
    xlim([-1.2 2.1])
    drawnow
    hold off
    title('Pyramidal')
end

clear dem

% Cerebellar
%--------------------------------------------------------------------------
DEM.M(1).V = exp(4);        % Restore sensory precision
DEM.M(1).E.s = 3/4;         % Overestimate smoothness

dem = spm_ADEM(DEM);

spm_figure('GetWin','Figure 1'); clf
for i = 1:size(dem.pU.x{1},2)
    MDP_DEM_MD_plot_arm(full(dem.pU.x{1}(1:3,i)),dem.M(1).pE.b1,dem.M(1).pE.b2,[-1;-1;-1]), hold on
    light
    view(10,0)
    zlim([-4 1])
    xlim([-1.2 2.1])
    drawnow
    hold off
    title('Cerebellar')
end

clear all

% COORDINATION
%==========================================================================
% Demonstration 2 - Mixed models (healthy)


% Specify discrete generatve model and mapping between continuous and
% discrete
%--------------------------------------------------------------------------

% Locations
%--------------------------------------------------------------------------

% 3D coordinates for three targets
L1  =   [   1;   -1; -3];
L2  =   [  -1;    1; -3];
L3  =   [   0;    0; -1];

% Intermediate (fictive) attracting points)
for i = 1:4
    L{i} = L1 + (i-1)*(L2 - L1)/4;
end
for i = 5:8
    L{i} = L2 + (i-5)*(L3 - L2)/4;
end
for i = 9:12
    L{i} = L3 + (i-9)*(L1 - L3)/4;
end


% MDP outcomes to DEM causes
%--------------------------------------------------------------------------
N  = 14;
nl = length(L);
nh = 3;

for i = 1:nh
    for j = 1:nl
        for k = 1:2
            c           = [L{j}; sparse(i,1,1,nh,1)];
            u           = [L{j}; sparse(i,1,1,nh,1)];
            demi.U{i,j,k} = u*ones(1,N);
            demi.C{i,j,k} = c*ones(1,N);
        end
    end
end

o     = [1 1 1];
O{1}  = spm_softmax(ones(nh,1));
O{2}  = spm_softmax(ones(nl,1));
O{3}  = spm_softmax(ones(2,1));

% generative model (Continuous)
%==========================================================================
DEM      = DEM_MDP_Movement_Disorders_GM;
DEM      = spm_MDP_DEM(DEM,demi,O,o);

% generative model (Discrete)
%==========================================================================
mdp      = MDP_DEM_Movement_Disorders_GM_1;
mdp.demi = demi;
mdp.DEM  = DEM;

mdp      = MDP_DEM_Movement_Disorders_GM_2(mdp);

% invert or solve
%--------------------------------------------------------------------------
rng default
MDP  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(MDP);

% illustrate belief updating - discrete
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
spm_MDP_VB_LFP(MDP,[],2);

% illustrate belief updating - continuous (last movement)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf;
spm_DEM_qU(MDP.mdp(end).dem(end).qU)

% illustrate movement
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf

[Q, R, X, Y, Z] = MDP_DEM_Movement_Disorders_Animate(MDP,L);


spm_figure('GetWin','Figure 4'); clf
plot3(Q(1,:),Q(2,:),Q(3,:),'LineWidth',2,'Color','k'), hold on
for m = 1:length(MDP.mdp)
    for j = 1:length(MDP.mdp(m).dem)
        h = (m - 1)*length(MDP.mdp(m).dem) + j;
        H = length(MDP.mdp)*length(MDP.mdp(m).dem);
        surf(0.1*X + MDP.mdp(m).dem(j).qU.v{2}(1,i),0.1*Y + MDP.mdp(m).dem(j).qU.v{2}(2,i),0.1*Z + MDP.mdp(m).dem(j).qU.v{2}(3,2),'Facecolor','r','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'FaceAlpha',h/(1.5*H))
    end
end

light
axis equal
spm_figure('GetWin','Figure 5'); clf
subplot(5,1,1)
plot(0.05*(0:size(R,2)-1),R' - repmat(R(:,1)',size(R,2),1))
title('Normal')

% Demonstration 3 - Mixed models (Pyramidal versus extrapyramidal)
%--------------------------------------------------------------------------
% Pyramidal

mdp.MDP.DEM.M(1).V = exp(5); % Overestimate precision

% invert or solve
%--------------------------------------------------------------------------
rng default
MDP  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(MDP);

spm_figure('GetWin','Figure 6'); clf

[Q, R, X, Y, Z] = MDP_DEM_Movement_Disorders_Animate(MDP,L);


spm_figure('GetWin','Figure 6'); clf
plot3(Q(1,:),Q(2,:),Q(3,:),'LineWidth',2,'Color','k'), hold on
for m = 1:length(MDP.mdp)
    for j = 1:length(MDP.mdp(m).dem)
        h = (m - 1)*length(MDP.mdp(m).dem) + j;
        H = length(MDP.mdp)*length(MDP.mdp(m).dem);
        surf(0.1*X + MDP.mdp(m).dem(j).qU.v{2}(1,i),0.1*Y + MDP.mdp(m).dem(j).qU.v{2}(2,i),0.1*Z + MDP.mdp(m).dem(j).qU.v{2}(3,2),'Facecolor','r','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'FaceAlpha',h/(1.5*H))
    end
end
light
axis equal

spm_figure('GetWin','Figure 5');
subplot(5,1,2)
plot(0.05*(0:size(R,2)-1),R' - repmat(R(:,1)',size(R,2),1))
title('Pyramidal')

% Extrapyramidal
%--------------------------------------------------------------------------
mdp.MDP.DEM.M(1).V = exp(4);
mdp.MDP.beta = 50;

% invert or solve
%--------------------------------------------------------------------------
rng default
MDP  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(MDP);

% illustrate movement
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 8'); clf

[Q, R, X, Y, Z] = MDP_DEM_Movement_Disorders_Animate(MDP,L);


spm_figure('GetWin','Figure 8'); clf
plot3(Q(1,:),Q(2,:),Q(3,:),'LineWidth',2,'Color','k'), hold on
for m = 1:length(MDP.mdp)
    for j = 1:length(MDP.mdp(m).dem)
        h = (m - 1)*length(MDP.mdp(m).dem) + j;
        H = length(MDP.mdp)*length(MDP.mdp(m).dem);
        surf(0.1*X + MDP.mdp(m).dem(j).qU.v{2}(1,i),0.1*Y + MDP.mdp(m).dem(j).qU.v{2}(2,i),0.1*Z + MDP.mdp(m).dem(j).qU.v{2}(3,2),'Facecolor','r','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'FaceAlpha',h/(1.5*H))
    end
end
light
axis equal

spm_figure('GetWin','Figure 5');
subplot(5,1,3)
plot(0.05*(0:size(R,2)-1),R' - repmat(R(:,1)',size(R,2),1))
title('Extrapyramidal')

% Cerebellar
%--------------------------------------------------------------------------
mdp.MDP.beta = 16;
mdp.MDP.DEM.M(1).E.s = 3/4;

% invert or solve
%--------------------------------------------------------------------------
rng default
MDP  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(MDP);

% illustrate movement
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 9'); clf

[Q, R, X, Y, Z] = MDP_DEM_Movement_Disorders_Animate(MDP,L);

spm_figure('GetWin','Figure 9'); clf
plot3(Q(1,:),Q(2,:),Q(3,:),'LineWidth',2,'Color','k'), hold on
for m = 1:length(MDP.mdp)
    for j = 1:length(MDP.mdp(m).dem)
        h = (m - 1)*length(MDP.mdp(m).dem) + j;
        H = length(MDP.mdp)*length(MDP.mdp(m).dem);
        surf(0.1*X + MDP.mdp(m).dem(j).qU.v{2}(1,i),0.1*Y + MDP.mdp(m).dem(j).qU.v{2}(2,i),0.1*Z + MDP.mdp(m).dem(j).qU.v{2}(3,2),'Facecolor','r','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'FaceAlpha',h/(1.5*H))
    end
end
light
axis equal

spm_figure('GetWin','Figure 5');
subplot(5,1,4)
plot(0.05*(0:size(R,2)-1),R' - repmat(R(:,1)',size(R,2),1))
title('Cerebellar')

% Executive
%--------------------------------------------------------------------------
mdp.MDP.DEM.M(1).E.s = 1/2;
mdp.A{3} = ones(size(mdp.A{3})); % Disconnect higher level

% invert or solve
%--------------------------------------------------------------------------
rng default
MDP  = spm_MDP_check(mdp);
MDP  = spm_MDP_VB_X(MDP);

% illustrate movement
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 10'); clf

[Q, R, X, Y, Z] = MDP_DEM_Movement_Disorders_Animate(MDP,L);

spm_figure('GetWin','Figure 10'); clf
plot3(Q(1,:),Q(2,:),Q(3,:),'LineWidth',2,'Color','k'), hold on
for m = 1:length(MDP.mdp)
    for j = 1:length(MDP.mdp(m).dem)
        h = (m - 1)*length(MDP.mdp(m).dem) + j;
        H = length(MDP.mdp)*length(MDP.mdp(m).dem);
        surf(0.1*X + MDP.mdp(m).dem(j).qU.v{2}(1,i),0.1*Y + MDP.mdp(m).dem(j).qU.v{2}(2,i),0.1*Z + MDP.mdp(m).dem(j).qU.v{2}(3,2),'Facecolor','r','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'FaceAlpha',h/(1.5*H))
    end
end
light
axis equal

spm_figure('GetWin','Figure 5');
subplot(5,1,5)
plot(0.05*(0:size(R,2)-1),R' - repmat(R(:,1)',size(R,2),1))
title('Executive')

function MDP = MDP_DEM_Movement_Disorders_GM_1

% First level
%--------------------------------------------------------------------------

% Initial states
%--------------------------------------------------------------------------
D{1} = ones(3,1); % Location of visual cue
D{2} = zeros(12,1);  % Location of arm
D{2}(1) = 1;

% Likelihood
%--------------------------------------------------------------------------
for f1 = 1:numel(D{1})
    for f2 = 1:numel(D{2})
        A{1}(f1,f1,f2) = 1; % Vision (target)
        A{2}(f2,f1,f2) = 1; % Vision + Proprioception (arm)
        A{3}(2,f1,f2)  = (f1 == ((f2-1)/4)+1);  % At target
        A{3}(1,f1,f2)  = ~(f1 == ((f2-1)/4)+1); % Not at target
    end
end

% Transitions
%--------------------------------------------------------------------------
B{1} = eye(3);

for k = 1:length(D{2})
    for j = 1:length(D{2})
        if k == j
            B{2}(k,j,1) = 1; % Stay still
        elseif (k - j) == 1
            B{2}(k,j,2) = 1; % Move to next point
        elseif (j - k) == 1
            B{2}(k,j,3) = 1; % Move to previous point
        end
    end
end

B{2}(1,end,2) = 1;
B{2}(end,1,3) = 1;

% Preferences
%--------------------------------------------------------------------------
C{1} = zeros(3,1);
C{2} = zeros(12,1);
C{3} = [-3 3]';

% Policies
%--------------------------------------------------------------------------
V(:,:,1) = ones(8,3); % (time, policy, factor)
V(:,:,2) = [1 2 3;
    1 2 3;
    1 2 3;
    1 2 3;
    1 1 1;
    1 1 1;
    1 1 1;
    1 1 1];

E = [1 0.5 0.5]';

% First level MDP
%--------------------------------------------------------------------------
MDP.A = A;
MDP.B = B;
MDP.C = C;
MDP.D = D;
MDP.E = E;
MDP.V = V;
MDP.T = 9;
MDP.chi = -exp(64);
MDP.beta = 16;
% MDP.beta = 50;
MDP     = spm_MDP_check(MDP);

function MDP = MDP_DEM_Movement_Disorders_GM_2(mdp)


% First level
%--------------------------------------------------------------------------
MDP.MDP = mdp;

% Second level
%--------------------------------------------------------------------------

% Initial states
%--------------------------------------------------------------------------
D{1} = ones(3,1); % Target location
D{2} = zeros(9,1); % Initial locations x policies at lower level ({L1,P1},{L1,P2},{L1,P3},{L2,P1},...)
D{2}([1 2 3]) = 1;

% Likelihood
%--------------------------------------------------------------------------
for f1 = 1:length(D{1})
    for f2 = 1:length(D{2})
        A{1}(f1,f1,f2) = 1; % Target location 1st level
        
        if f2 < 4
            A{2}(11,f1,f2) = 0.03;  % Start location 1st level
            A{2}(12,f1,f2) = 0.12;
            A{2}(1,f1,f2)  = 0.7;
            A{2}(2,f1,f2)  = 0.12;
            A{2}(3,f1,f2)  = 0.03;
            
        elseif f2 <7
            A{2}(3,f1,f2) = 0.03;
            A{2}(4,f1,f2) = 0.12;
            A{2}(5,f1,f2)  = 0.7;
            A{2}(6,f1,f2)  = 0.12;
            A{2}(7,f1,f2)  = 0.03;
        else
            A{2}(7,f1,f2) = 0.03;
            A{2}(8,f1,f2) = 0.12;
            A{2}(9,f1,f2)  = 0.7;
            A{2}(10,f1,f2)  = 0.12;
            A{2}(11,f1,f2)  = 0.03;
        end
        
        A{3}(1,f1,f2)  = (rem(f2+2,3) == 0) + 0.5;  % Policy 1st level (slight bias towards 'stay still' policy)
        A{3}(2,f1,f2)  = (rem(f2+1,3) == 0);
        A{3}(3,f1,f2)  = (rem(f2,3) == 0);
    end
end


% Transitions
%--------------------------------------------------------------------------
B{1}(:,:,1) = ones(length(D{1}));
B{1}(:,:,2) = ones(length(D{1})); % Action that causes no change to prevent entering HMM mode
I = eye(3);
B{2} = [I, [I(3,:);I(1:2,:)],([I(3,:);I(1:2,:)])']/3;
B{2} = kron(B{2},ones(3,1));

% Preferences
%--------------------------------------------------------------------------
C{1} = zeros(3,1);
C{2} = zeros(12,1);
C{3} = zeros(3,1);

% Second level MDP
%--------------------------------------------------------------------------
MDP.A = A;
MDP.B = B;
MDP.C = C;
MDP.D = D;
MDP.T = 4;

MDP.link  = [1 0 0;
    0 1 0];
MDP.linkE = [0 0 1];

MDP     = spm_MDP_check(MDP);

function DEM = DEM_MDP_Movement_Disorders_GM
% hidden states
%--------------------------------------------------------------------------
x   = [0;-1;0.5;0;0;0]; % shoulder rotation, flexion, elbow flexion (position and velocity)

% and actions
%--------------------------------------------------------------------------
a   = [0;0;0];

% parameters
%--------------------------------------------------------------------------
P.kv = 2;              % damping
P.a = a;
P.b1 = 2;
P.b2 = 1.5;

M(1).pE = P;

% recognition model
%==========================================================================
M(1).E.s = 1/2;                               % smoothness: default 1/2 [3/4 for lesion]
M(1).E.n = 4;                                 % order of
M(1).E.d = 2;                                 % generalised motion

% level 1
%--------------------------------------------------------------------------
M(1).f  = @(x,v,P) Mf1(x,v,P);                % plant dynamics
M(1).g  = @(x,v,P) Mg1(x,v,P);                % prediction

M(1).x  = x;                                  % hidden states
M(1).V  = exp(4);                             % error precision (g)[4] [3/8 ?]
M(1).W  = exp(2);                             % error precision (f) [2]

% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0;0;0;0;0;0];                      % priors (attracting point)
M(2).V  = exp(15);


% generative model
%==========================================================================

G(1).E.s = 1/2;                               % smoothness: default 1/2
G(1).E.n = 4;                                 % order of
G(1).E.d = 2;                                 % generalised motion

% first level
%--------------------------------------------------------------------------
G(1).f  = @(x,v,a,P) Gf1(x,v,a,P);
G(1).g  = @(x,v,a,P) Gg1(x,v,a,P);
G(1).x  = x;                                  % hidden states
G(1).V = exp(16);
% G(1).V  = exp([4*ones(1,6) 16*ones(1,6)]);    % error precision [16]
G(1).W  = exp(16);                            % error precision
G(1).pE = P;

% second level
%--------------------------------------------------------------------------
G(2).v  = [0;0;0;0;0;0];                      % exogenous forces
G(2).a  = spm_vec(a);                         % action forces
G(2).V  = exp(16);


DEM.G = G;
DEM.M = M;

% Equations for DEM model and process
%==========================================================================
function f = Mf1(x,v,P)
a = P.b1; % length of upper arm
b = P.b2; % length of forearm
c = [-1;-1;-1]; % location of shoulder

f = [x(4:6);MDP_DEM_MD_dthetadt(v(1:3),x,a,b,c) - P.kv*x(4:6)];

function g = Mg1(x,v,P)
a = P.b1; % length of upper arm
b = P.b2; % length of forearm
c = [-1;-1;-1]; % location of shoulder

g = [x;MDP_DEM_MD_transform(x(1:3),a,b,c);v(4:6)]; % proprioceptive, visual (hand), visual (target)

function f = Gf1(x,~,a,P)
f = [x(4:6);a  - P.kv*x(4:6)];

function g = Gg1(x,v,~,P)
a = P.b1; % length of upper arm
b = P.b2; % length of forearm
c = [-1;-1;-1]; % location of shoulder

g = [x;MDP_DEM_MD_transform(x(1:3),a,b,c);v(4:6)]; % proprioceptive, visual (hand), visual (target)

function g = Gg1Reflex(x,v,~,P)
a = P.b1; % length of upper arm
b = P.b2; % length of forearm
c = [-1;-1;-1]; % location of shoulder

g = [x + [zeros(3,1);v(1:3)];MDP_DEM_MD_transform(x(1:3),a,b,c);v(4:6)]; % proprioceptive, visual (hand), visual (target)

% Additional routines (plotting/coordinate transforms)
%==========================================================================
function y = MDP_DEM_MD_transform(x,a,b,c)
% converts angles to hand position in euclidean space
% x = angle of shoulder rotation, flexion, elbow flexion

y1 = c + [a*cos(x(2))*cos(x(1));
    a*cos(x(2))*sin(x(1));
    a*sin(x(2))];

y = y1 + [b*cos(x(3)+x(2))*cos(x(1));
    b*cos(x(3)+x(2))*sin(x(1));
    b*sin(x(3)+x(2))];

function MDP_DEM_MD_plot_arm(x,a,b,c)

y1 = c + [a*cos(x(2))*cos(x(1));
    a*cos(x(2))*sin(x(1));
    a*sin(x(2))];

y = y1 + [b*cos(x(3)+x(2))*cos(x(1));
    b*cos(x(3)+x(2))*sin(x(1));
    b*sin(x(3)+x(2))];

[S1, S2, S3] = sphere;

surf(S1*0.2+c(1),S2*0.2+c(2),S3*0.2+c(3),'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7), hold on
surf(S1*0.1+y1(1),S2*0.1+y1(2),S3*0.1+y1(3),'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7)
surf(S1*0.07+y(1),S2*0.07+y(2),S3*0.07+y(3),'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7)

plot3([c(1) y1(1)],[c(2) y1(2)],[c(3) y1(3)],'LineWidth',6,'Color',[0.1 0.1 0.8])
plot3([y1(1) y(1)],[y1(2) y(2)],[y1(3) y(3)],'LineWidth',5,'Color',[0.1 0.1 0.8])
plot3([y1(1) y(1)],[y1(2) y(2)],[y1(3)-0.1 y(3)-0.1],'LineWidth',4,'Color',[0.1 0.1 0.8])

MDP_DEM_Movements_Hand([pi,pi/2 + atan((y(3)-y1(3))/(sqrt((y(1) - y1(1))^2+(y(2) - y1(2))^2))),x(1)],0.25,y)
axis equal

function y = MDP_DEM_MD_dthetadt(v,x,a,b,c)

y1 = c + [a*cos(x(2))*cos(x(1));
    a*cos(x(2))*sin(x(1));
    a*sin(x(2))];

y2 = y1 + [b*cos(x(3)+x(2))*cos(x(1));
    b*cos(x(3)+x(2))*sin(x(1));
    b*sin(x(3)+x(2))];


y(1) = [v - y2]'*[-((a*cos(x(2)) + b*cos(x(3)+x(2))*sin(x(1))));
    ((a*cos(x(2)) + b*cos(x(3)+x(2))*cos(x(1))));
    0];

y(2) = [v - y2]'*[-((a*sin(x(2)) + b*sin(x(3) + x(2)))*cos(x(1)));
    -((a*sin(x(2)) + b*sin(x(3) + x(2))*sin(x(1))));
    (a*cos(x(2)) +(b*cos(x(2) + x(3))))];

y(3) = [v - y2]'*[-(a*sin(x(3) + x(2))*cos(x(1)));
    -(a*sin(x(3)+ x(2))*sin(x(1)));
    (b*cos(x(2) + x(3)))];


y = 0.05*y'; % 0.05

function MDP_DEM_Movements_Hand(a,s,l)
% angle    - a
% scale    - s
% location - l

% Carpals (as a single sphere)
%--------------------------------------------------------------------------
[X,Y,Z] = sphere;
Rx = makehgtform('xrotate',a(1));
Ry = makehgtform('yrotate',a(2));
Rz = makehgtform('zrotate',a(3));
T  = makehgtform('translate',l);
S  = makehgtform('scale',s);
q = hgtransform;

LW = 5;

% 2nd finger
%--------------------------------------------------------------------------

% MCP
%--------------------------------------------------------------------------
plot3([0 0],[0 0],[0 1.1],'b','LineWidth',LW,'Parent',q)
surf(0.1*X, 0.1*Y ,0.1*Z + 1.1,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% PIP
%--------------------------------------------------------------------------
plot3([0 0],[0 0.4],[1.1 1.4],'b','LineWidth',LW,'Parent',q)
surf(0.1*X, 0.1*Y + 0.4 ,0.1*Z + 1.4,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% MIP
%--------------------------------------------------------------------------
plot3([0 0],[0.4 0.7],[1.4 1.5],'b','LineWidth',LW,'Parent',q)
surf(0.1*X, 0.1*Y + 0.7 ,0.1*Z + 1.5,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% DIP
%--------------------------------------------------------------------------
plot3([0 0],[0.7 1],[1.5 1.6],'b','LineWidth',LW,'Parent',q)
surf(0.1*X, 0.1*Y + 1 ,0.1*Z + 1.6,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)


% 3rd finger
%--------------------------------------------------------------------------

% MCP
%--------------------------------------------------------------------------
plot3([-0.1 -0.2],[0 0],[0 1.05],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.2, 0.1*Y ,0.1*Z + 1.05,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% PIP
%--------------------------------------------------------------------------
plot3([-0.2 -0.2],[0 0.3],[1.05 1.4],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.2, 0.1*Y + 0.3 ,0.1*Z + 1.4,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% MIP
%--------------------------------------------------------------------------
plot3([-0.2 -0.2],[0.3 0.6],[1.4 1.5],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.2, 0.1*Y + 0.6 ,0.1*Z + 1.5,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% DIP
%--------------------------------------------------------------------------
plot3([-0.2 -0.2],[0.6 0.9],[1.5 1.6],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.2, 0.1*Y + 0.9 ,0.1*Z + 1.6,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)


% 4th finger
%--------------------------------------------------------------------------

% MCP
%--------------------------------------------------------------------------
plot3([-0.2 -0.4],[0 0],[0 1],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.4, 0.1*Y ,0.1*Z + 1,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% PIP
%--------------------------------------------------------------------------
plot3([-0.4 -0.4],[0 0.2],[1 1.4],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.4, 0.1*Y + 0.2 ,0.1*Z + 1.4,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% MIP
%--------------------------------------------------------------------------
plot3([-0.4 -0.4],[0.2 0.5],[1.4 1.5],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.4, 0.1*Y + 0.5 ,0.1*Z + 1.5,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% DIP
%--------------------------------------------------------------------------
plot3([-0.4 -0.4],[0.5 0.8],[1.5 1.6],'b','LineWidth',LW,'Parent',q)
surf(0.1*X - 0.4, 0.1*Y + 0.8 ,0.1*Z + 1.6,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)

% 1st finger
%--------------------------------------------------------------------------

% MCP
%--------------------------------------------------------------------------
plot3([0.1 0.2],[0 0],[0 1],'b','LineWidth',LW,'Parent',q)
surf(0.1*X + 0.2, 0.1*Y ,0.1*Z + 1,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% PIP
%--------------------------------------------------------------------------
plot3([0.2 0.2],[0 0.3],[1 1.4],'b','LineWidth',LW,'Parent',q)
surf(0.1*X + 0.2, 0.1*Y + 0.3 ,0.1*Z + 1.4,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% MIP
%--------------------------------------------------------------------------
plot3([0.2 0.2],[0.3 0.6],[1.4 1.5],'b','LineWidth',LW,'Parent',q)
surf(0.1*X + 0.2, 0.1*Y + 0.6 ,0.1*Z + 1.5,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% DIP
%--------------------------------------------------------------------------
plot3([0.2 0.2],[0.6 0.9],[1.5 1.6],'b','LineWidth',LW,'Parent',q)
surf(0.1*X + 0.2, 0.1*Y + 0.9 ,0.1*Z + 1.6,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)

% Thumb
%--------------------------------------------------------------------------

% MCP
%--------------------------------------------------------------------------
plot3([0.2 0.5],[0 0.05],[0 0.4],'b','LineWidth',LW,'Parent',q)
surf(0.1*X + 0.5, 0.1*Y +0.05,0.1*Z + 0.4,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% PIP
%--------------------------------------------------------------------------
plot3([0.5 0.5],[0.05 0.45],[0.4 0.6],'b','LineWidth',LW,'Parent',q)
surf(0.1*X + 0.5, 0.1*Y +0.45,0.1*Z + 0.6,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)
% DIP
%--------------------------------------------------------------------------
plot3([0.5 0.45],[0.45 0.6],[0.6 0.62],'b','LineWidth',LW,'Parent',q)
surf(0.1*X + 0.45, 0.1*Y +0.6,0.1*Z + 0.62,'Facecolor','b','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7,'Parent',q)

set(q,'Matrix',T*Rz*Rx*Ry*S);

function [Q, R, X, Y, Z] = MDP_DEM_Movement_Disorders_Animate(MDP,L)
opengl('software')
Q = [];
R = [];

for m = 1:length(MDP.mdp)
    for j = 1:length(MDP.mdp(m).dem)
        for i = 2:size(MDP.mdp(m).dem(j).pU.x{1},2)
            MDP_DEM_MD_plot_arm(full(MDP.mdp(m).dem(j).pU.x{1}(1:3,i)),MDP.mdp(m).dem(j).M(1).pE.b1,MDP.mdp(m).dem(j).M(1).pE.b2,[-1;-1;-1]), hold on
            [X,Y,Z] = sphere;
            surf(0.1*X + MDP.mdp(m).dem(j).qU.v{2}(1,i),0.1*Y + MDP.mdp(m).dem(j).qU.v{2}(2,i),0.1*Z + MDP.mdp(m).dem(j).qU.v{2}(3,i),'Facecolor','r','EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7)
            for k = 1:3
                surf(0.2*X + L{(k-1)*4+1}(1),0.2*Y + L{(k-1)*4+1}(2),0.2*Z + L{(k-1)*4+1}(3),'Facecolor',[1 1 1]-MDP.mdp(m).dem(j).Y(9+k,i),'EdgeColor','none','FaceLighting','gouraud','AmbientStrength',0.7)
            end
            axis equal
            light
            drawnow
            hold off
                       
            Q(:,end+1) = MDP_DEM_MD_transform(full(MDP.mdp(m).dem(j).pU.x{1}(1:3,i)),MDP.mdp(m).dem(j).M(1).pE.b1,MDP.mdp(m).dem(j).M(1).pE.b2,[-1;-1;-1]);
        end
        R = [R full(MDP.mdp(m).dem(j).pU.x{1}(1:3,:))];
    end
end