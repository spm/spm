% Action-DEM demo specifying an attractor (in terms of the parameters of
% desired equations of motion) This demo first exposes an agent to a Lorenz
% attractor world such that it can learn the three parameters of the Lorenz
% system. It is then placed in a test world with a fixed point attractor to
% see if it has remembered the chaotic dynamics it learnt in the training
% environment (in DEMO mode previously optimised models are loaded)


% generative model
%==========================================================================
clear
DEMO     = 1;                          % switch for demo
G(1).E.s = 1/4;                        % smoothness
G(1).E.n = 6;                          % smoothness
G(1).E.d = 2;                          % smoothness
G(1).E.linear = 3;                     % nonlinear


% dynamics
%--------------------------------------------------------------------------
f0      = '[v; 0; 0] + [-1 1 1; -1/2 -1 0; 0 1 -1]*x/P';
fA      = '[v; 0; 0] + [-1 1 1; -1/2 -1 0; 0 1 -1]*x/P + a';
fL      = '[v; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/64';

% parameters
%--------------------------------------------------------------------------
P0      = 16;
PL      = [10; -8/3; 32];
PM      = [0;     0;  0];

% P(1): Prandtl number
% P(2): 8/3
% P(3): Rayleigh number

% level 1
%--------------------------------------------------------------------------
G(1).x  = [1; 1; 24];
G(1).f  = inline(fL ,'x','v','a','P');
G(1).g  = inline('x','x','v','a','P');
G(1).pE = PL;
G(1).V  = exp(8);                           % error precision
G(1).W  = exp(8);                           % error precision

% level 2
%--------------------------------------------------------------------------
G(2).a  = [0;0;0];                          % action variables
G(2).v  = 0;                                % inputs
G(2).V  = exp(16);
G       = spm_ADEM_M_set(G);

% plot flow fields and equilibrium densities
%==========================================================================
spm_figure('GetWin','Graphics');

x{1}    = linspace(-20,20,32);
x{2}    = linspace(-32,32,32);
x{3}    = linspace(  10,40,8);


% controlled flow (P0)
%--------------------------------------------------------------------------
subplot(3,2,1)
spm_fp_display_density(G,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('controlled','Fontsize',16)
a     = axis;

% exemplar trajectory
%--------------------------------------------------------------------------
U.u   = sparse(1024,G(1).m);
t     = spm_int_J(G(1).pE,G,U);
subplot(3,2,2)
plot(t(:,1),t(:,2),'r')
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('controlled','Fontsize',16)
axis(a);


% uncontrolled flow (P0)
%--------------------------------------------------------------------------
G0       = G;
G0(1).f  = inline(f0,'x','v','P');
G0(1).pE = P0;

subplot(3,2,3)
spm_fp_display_density(G0,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('before','Fontsize',16)
a = axis;

% exemplar trajectory
%--------------------------------------------------------------------------
t     = spm_int_J(G0(1).pE,G0,U);
subplot(3,2,4)
plot(t(:,1),t(:,2),'r')
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('before','Fontsize',16)
axis(a);
drawnow


% recognition model: learn the controlled environmental dynamics
%==========================================================================


% make a naive model (M)
%--------------------------------------------------------------------------
M       = G;
M(1).f  = inline(fL ,'x','v','P');
M(1).g  = inline('x','x','v','P');
M(1).pE = PM;
M(1).pC = eye(3,3)*4;


% teach naive model by exposing it to a controlled environment (G)
%--------------------------------------------------------------------------
clear DEM

% length of realization
%--------------------------------------------------------------------------
n     = 256;
C     = sparse(n,1);

DEM.M = M;
DEM.G = G;
DEM.C = C;
DEM.U = C;

% optimise recognition model
%--------------------------------------------------------------------------
if DEMO
    for i = 1:16
        DEM = spm_ADEM(DEM);
        DEM.M(1).pE = DEM.qP.P{1};
        DEM.M(1).x  = DEM.qU.x{1}(:,end);
        DEM.M(2).a  = DEM.qU.a{2}(:,end);
    end
    
    save DEM_lorenz G DEM
else
    try
        load DEM_lorenz
    catch
        disp({'sorry; the previously prepared .mat file is not compatible';
              'with your version of matlab - run this routine with DEMO = 1'})
        return
    end
end

spm_figure('GetWin','DEM');
spm_DEM_qP(DEM.qP,DEM.pP)
spm_DEM_qU(DEM.qU)


% replace priors with learned conditional expectation and plot
%--------------------------------------------------------------------------
M             = DEM.M;
M(1).E.linear = 3;

spm_figure('GetWin','Graphics');
subplot(3,2,5)
spm_fp_display_density(M,x);
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('after','Fontsize',16)
a   = axis;

% exemplar trajectory
%--------------------------------------------------------------------------
t   = spm_int_J(M(1).pE,M,U);
subplot(3,2,6)
plot(t(:,1),t(:,2),'r')
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('after','Fontsize',16)
axis(a);
drawnow


% evaluate performance under active inference
%==========================================================================
N       = 256;
C       = spm_conv(randn(1,N)*4,8);

% prescribe desired motion directly with true parameters
%--------------------------------------------------------------------------
M(1).pE = PL;
M(1).pC = sparse(3,3);

% create uncontrolled environment: A (with action)
%--------------------------------------------------------------------------
A       = G;
A(1).f  = inline(f0 ,'x','v','a','P');       % no action
A(1).f  = inline(fA ,'x','v','a','P');       % with action
A(1).g  = inline('x','x','v','a','P');
A(1).pE = P0;
A(1).pC = 0;

% make the recognition model confident about its predictions
%--------------------------------------------------------------------------
M(1).x  = G(1).x;
M(1).W  = exp(12);
M(1).V  = exp(8);
A(1).W  = exp(16);

% create DEM structure
%--------------------------------------------------------------------------
clear DEM
DEM.G   = A;
DEM.M   = M;
DEM.C   = sparse(1,N);                       % no perturbation
DEM.C   = C;                                 % with perturbation
DEM.U   = sparse(1,N);
DEM     = spm_ADEM(DEM);

spm_DEM_qU(DEM.qU,DEM.pU)

% plot behaviour
%--------------------------------------------------------------------------
subplot(2,2,3)
spm_fp_display_density(M,x); hold on
plot(DEM.pU.v{1}(1,:),DEM.pU.v{1}(2,:),'b'), hold off
xlabel('position','Fontsize',12)
ylabel('velocity','Fontsize',12)
title('learnt','Fontsize',16)





return

% Notes: alternative specification of Lorenz system
%==========================================================================

% dynamics
%--------------------------------------------------------------------------
fL      = '[v; 0; 0] + [-P(1) P(1) 0; P(3) -1 -x(1); x(2) 0 P(2)]*x/32';

% parameters
%--------------------------------------------------------------------------
PL      = [10; -8/3; 32];
P.A     = [-PL(1) PL(1) 0; PL(3) -1 0; 0 0 PL(2)]/32;
P.B     = [0 0  0;
    0 0 -1;
    0 1  0]/32;

% PL(1): Prandtl number
% PL(2): 8/3
% PL(3): Rayleigh number
