function ADEM_occulomotor_delays
% Oculomotor following and delays under active inference:
%__________________________________________________________________________
% This demo illustrates oculomotor following and slow pursuit. The focus
% here is on oculomotor delays and their compensation in generalised
% coordinates of motion. This is illustrates using a ‘sweep’ of motion to
% examine the effects of motor delays, sensory delays and their interaction
% under active inference. We then move on to oculomotor following of sine
% wave motion, in which the trajectory entrains following under compensated
% dynamics. This entrainment can be destroyed by rectifying the sine wave
% which calls for a more realistic (hierarchical) model of motion that
% registers its phase and anticipates the next onset of motion (cf
% movement behind an occluder). These simulations depend delicately on
% delays and precisions (gains) that are chosen to make oculomotor following
% under uncertainty relatively slow. The dependency on visual uncertainty
% (contrast) is illustrated by changing the precision of the generalised
% motion of the sensory target.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_occulomotor_delays.m 4804 2012-07-26 13:14:18Z karl $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states:
%   x.o(1) - oculomotor angle
%   x.o(2) - oculomotor angle
%   x.x(1) - target location (visual) - extrinsic coordinates (Cartesian)
%   x.x(2) - target location (visual) - extrinsic coordinates (Cartesian)
%
% v    - causal states: force on target
%
% g    - sensations:
%   g(1) - oculomotor angle (proprioception)
%   g(2) - oculomotor angle (proprioception)
%   g(3) - target location (visual) - intrinsic coordinates (polar)
%   g(4) - target location (visual) - intrinsic coordinates (polar)
%--------------------------------------------------------------------------
 
 
% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
x.o = [0;0];                                  % oculomotor angle
x.x = [0;0];                                  % target location
 
% Recognition model
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 4;                                 % order of
M(1).E.d = 1;                                 % generalised motion
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = '[(x.x - x.o)/2; v - x.x]';        % extrinsic coordinates
M(1).g   = '[x.o; x.x - x.o]';                % intrinsic coordinate
M(1).x   = x;                                 % hidden states
M(1).V   = exp(4);                            % error precision
M(1).W   = exp(4);                            % error precision
 
% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0; 0];                             % inputs
M(2).V  = exp(0);
 
 
% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = '[a/4 - x.o/32; v - x.x]';          % extrinsic coordinates
G(1).g  = '[x.o; x.x - x.o]';                 % intrinsic coordinates
G(1).x  = x;                                  % hidden states
G(1).V  = exp(16);                            % error precision (errors)
G(1).W  = exp(16);                            % error precision (motion)
G(1).U  = [1 1 0 0]*exp(2);                   % error precision (action)
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0];                             % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);
 
 
% Solve or integrate
%==========================================================================
N       = 64;                                 % length of data sequence
dt      = 16;                                 % time step (ms)
t       = (1:N)*dt;                           % PST
C       = [1; 0]*spm_phi(((1:N) - N/3)/2);    % simple motion (sweep)
DEM.G   = G;
DEM.M   = M;
DEM.C   = C;
DEM.U   = zeros(2,N);
 
% compensated delays
%--------------------------------------------------------------------------
DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;
DEM         = spm_ADEM(DEM);
 
% uncompensated sensory delays
%--------------------------------------------------------------------------
DEM.M(1).Ta = 0;
DEM.M(1).Ty = 2;
DEMy        = spm_ADEM(DEM);
 
% uncompensated motor delays
%--------------------------------------------------------------------------
DEM.M(1).Ta = 2;
DEM.M(1).Ty = 0;
DEMa        = spm_ADEM(DEM);
 
% uncompensated delays
%--------------------------------------------------------------------------
DEM.M(1).Ta = 2;
DEM.M(1).Ty = 2;
DEMay       = spm_ADEM(DEM);
 
 
% Compare latencies graphically
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(2,2,1), title({'Oculomotor following model', 'prediction and error'},'FontSize',16)
subplot(2,2,2), title({'Simple (sweep) motion', 'hidden states'},'FontSize',16)
 
 
spm_figure('GetWin','Figure 2'); clf
subplot(3,2,1); plot(t,DEM.pU.v{1},'b',t,DEM.qU.v{1},':b'), hold on
subplot(3,2,2); plot(t,DEM.qU.a{2},'b'), hold on
subplot(3,2,1); plot(t,DEMy.pU.v{1},'r',t,DEMy.qU.v{1},':r'), hold on
subplot(3,2,2); plot(t,DEMy.qU.a{2},'r'), hold on
 
subplot(3,2,3); plot(t,DEM.pU.v{1},'b',t,DEM.qU.v{1},':b'), hold on
subplot(3,2,4); plot(t,DEM.qU.a{2},'b'), hold on
subplot(3,2,3); plot(t,DEMa.pU.v{1},'r',t,DEMa.qU.v{1},':r'), hold on
subplot(3,2,4); plot(t,DEMa.qU.a{2},'r'), hold on
 
subplot(3,2,5); plot(t,DEM.pU.v{1},'b',t,DEM.qU.v{1},':b'), hold on
subplot(3,2,6); plot(t,DEM.qU.a{2},'b'), hold on
subplot(3,2,5); plot(t,DEMay.pU.v{1},'r',t,DEMay.qU.v{1},':r'), hold on
subplot(3,2,6); plot(t,DEMay.qU.a{2},'r'), hold on
 
 
% labels and titles
%--------------------------------------------------------------------------
subplot(3,2,1); title('Sensory input and predictions','FontSize',16), spm_axis tight
subplot(3,2,2); title('Target and eye positions','FontSize',16), spm_axis tight
 
subplot(3,2,5); xlabel('time (ms)','FontSize',12), spm_axis tight
subplot(3,2,6); xlabel('time (ms)','FontSize',12), spm_axis tight
subplot(3,2,4); spm_axis tight
 
subplot(3,2,1); ylabel('Sensory delays','FontSize',16), spm_axis tight
subplot(3,2,3); ylabel('Motor delays','FontSize',16), spm_axis tight
subplot(3,2,5); ylabel('Sensorimotor delays','FontSize',16), spm_axis tight
 
 
% Contrast-dependent latencies (compensated for delays)
%==========================================================================
DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;
DEM.G(1).U  = [1 1 0 0]*exp(4);
 
% change prior precision (contrast)
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
 
V = (2:8);                         
for i = 1:length(V)
    DEM.M(1).V = diag(exp([4 4 V(i) V(i)]));
    DEM        = spm_ADEM(DEM);
    pd         = sqrt(sum(DEM.pU.v{1}(3:4,:).^2));
    qd         = sqrt(sum(DEM.qU.v{1}(3:4,:).^2));
    pD(i)      = max(pd);
    qD(i)      = max(qd);
    
    spm_figure('GetWin','Figure 3');
    subplot(2,1,1)
    plot(t,pd,'k',t,qd,':k'), hold on
    
end
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');
subplot(2,1,1), hold off
xlabel('Contrast','FontSize',12)
ylabel('Spatial lag','FontSize',12)
title('Distance from target','FontSize',16)
axis square
 
subplot(2,1,2)
plot(V,pD,'-k',V,pD,'ok',V,qD,':k',V,qD,'xk')
xlabel('Contrast','FontSize',12)
ylabel('Spatial lag','FontSize',12)
title('Peak lag','FontSize',16)
axis square
 
 
 
% Oculomotor latencies (sinusoidal movement)
%==========================================================================
DEM.M(1).Ta = 0;
DEM.M(1).Ty = 0;
DEM.G(1).U  = [1 1 0 0]*exp(8);
 
% oculomotor following with (first order) generative model and sine wave
%--------------------------------------------------------------------------
C       = [-1; 0]*(sin((1:N)*2*pi/32).*((1:N) > 16));
DEM.C   = C;
DEM.U   = zeros(2,N);
DEM     = spm_ADEM(DEM);
 
spm_figure('GetWin','Figure 4');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(2,2,1), title({'Oculomotor following:', 'prediction and error'},'FontSize',16)
subplot(2,2,2), title({'Periodic motion:', 'hidden states'},'FontSize',16)
 
 
% slow pursuit following with (second order) generative model
%--------------------------------------------------------------------------
mx.o = [0;0];                                  % oculomotor angle
mx.x = [0;0];                                  % target location
mx.v = [0;0];                                  % oculomotor velocity
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = '[x.v; v - x.x; (v - x.o)/4 - x.v/2]';
M(1).g   = '[x.o; x.x - x.o]';                 % intrinsic coordinate
M(1).x   = mx;                                 % hidden states
M(1).V   = exp(4);                             % error precision
M(1).W   = exp(4);                             % error precision
 
% level 2:
%--------------------------------------------------------------------------
M(2).v  = [0; 0];                              % inputs
M(2).V  = exp(0);
 
 
% generative model
%==========================================================================
gx.o = [0;0];                                  % oculomotor angle
gx.x = [0;0];                                  % target location
gx.v = [0;0];                                  % oculomotor velocity
 
% first level
%--------------------------------------------------------------------------
G(1).f  = '[x.v; v - x.x; a/4 - x.v]';
G(1).g  = '[x.o; x.x - x.o]';                 % intrinsic coordinates
G(1).x  = gx;                                 % hidden states
G(1).V  = exp(16);                            % error precision (errors)
G(1).W  = exp(16);                            % error precision (motion)
G(1).U  = [1 1 0 0]*exp(8);                   % motor gain
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0];                             % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);
 
% Sine wave cause
%--------------------------------------------------------------------------
DEM.M   = M;
DEM.G   = G;
DEM.C   = C;
DEM.U   = zeros(2,N);
DEM     = spm_ADEM(DEM);
 
spm_figure('GetWin','Figure 5');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(2,2,1), title({'Slow pursuit:', 'prediction and error'},'FontSize',16)
subplot(2,2,2), title({'Periodic motion:', 'hidden states'},'FontSize',16)
 
 
% rectified sinusoidal movement
%--------------------------------------------------------------------------
DEM.C      = exp((C - 1)*4);
DEM        = spm_ADEM(DEM);
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 6');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(2,2,1), title({'Slow pursuit:', 'prediction and error'},'FontSize',16)
subplot(2,2,2), title({'Aperiodic motion:', 'hidden states'},'FontSize',16)
 
 
 
% Slow pursuit with a hierarchical generative model
%==========================================================================
% Endow the model with internal dynamics (a simple oscillator) so that is
% recognises and remembers the trajectory to anticipate jumps in rectified
% sinusoidal motion.
 
 
% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = '[x.v; v - x.x; (v - x.o)/4 - x.v/2]';
M(1).g   = '[x.o; x.x - x.o]';                 % intrinsic coordinate
M(1).x   = mx;                                 % hidden states
M(1).V   = exp(4);                             % error precision
M(1).W   = exp(4);                             % error precision
 
% level 2: With hidden (memory) states
%--------------------------------------------------------------------------
M(2).f  = '[0 1; -1 0]*x*v'; 
M(2).g  = '[exp((x(1) - 1)*4); 0]';          
M(2).x  = [0; 0];                              % hidden states
M(2).V  = exp(4);                              % error precision
M(2).W  = exp(4);                              % error precision
 
% level 3: Encoding frequency of memory states (U)
%--------------------------------------------------------------------------
U       = 2*pi/32;
M(3).v  = U;
M(3).V  = exp(4);
 
% give the model prior beliefs about angular motion (but not its phase)
%--------------------------------------------------------------------------
DEM.G   = G;
DEM.M   = M;
DEM.C   = exp((C - 1)*4);
DEM.U   = zeros(1,N) + U;
DEM     = spm_ADEM(DEM);
 
% Graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(3,2,1), title({'Anticipatory model:', 'prediction and error'},'FontSize',16)
subplot(3,2,2), title({'Aperiodic motion:', 'hidden states'},'FontSize',16)
 
 
% create movie in extrinsic and intrinsic coordinates
%==========================================================================
spm_figure('GetWin','Figure 8'); clf
spm_dem_pursuit_movie(DEM,0)
