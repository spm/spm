function ADEM_observe
% This demo illustrates action-observation using synthetic writing under 
% active inference. It shows how expectations about hidden states can be 
% both cause and consequence of observed action (of self and others 
% respectively). We first illustrate the generation of behaviour using a 
% Lotka-Volterra form stable heteroclinic orbit. We then reproduce the 
% same forces on the agent's arm but switching off the precision of 
% proprioceptive inputs. This can be seen as attending selectively to 
% visual inputs. The resulting inference calls upon the same hidden-states 
% and implicit predictions (in a generalised or dynamic sense). These 
% simulations can be regarded as simulations of mirror neuron responses.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_observe.m 7679 2019-10-24 15:54:07Z spm $
 
 
% hidden causes and states
%==========================================================================
% x    - hidden states
%   x.x(1) - joint angle
%   x.x(2) - joint angle
%   x.x(3) - angular velocity
%   x.x(4) - angular velocity
%
%   x.a(1) - attraction (location 1)
%   x.a(2) - attraction (location 2)
%   x.a(3) - attraction (location 3)
%    ...
%
% v    - causal states
%   v(1) - not used
%
%--------------------------------------------------------------------------


% parameters mapping from (unstable) point attractors to visual space
%--------------------------------------------------------------------------
P = [1.0  1.0;
     1.1  1.2;
     1.0  0.4;
     1.0  1.0;
     1.4  0.9;
     0.9  1.0]';
n = size(P,2);                                % number of attractors
 
 
    
% Recognition model (linear for expediency)
%==========================================================================
M(1).E.s = 1/2;                               % smoothness
M(1).E.n = 4;                                 % order of 
M(1).E.d = 2;                                 % generalised motion


% level 1: Displacement dynamics and mapping to sensory/proprioception
%--------------------------------------------------------------------------
M(1).f   = 'spm_fx_dem_observe';              % plant dynamics
M(1).g   = 'spm_gx_dem_write';                % prediction
 
M(1).x.x = [pi/2; pi/2; 0; 0];                % physical states
M(1).x.a = sparse(1,1,3,n,1) - 4;             % attractor states
M(1).V   = exp(4);                            % error precision
M(1).W   = exp(8);                            % error precision
M(1).pE  = P;
 
 
% level 2: not used
%--------------------------------------------------------------------------
M(2).v  = 0;                                  % inputs
M(2).V  = exp(8);
 
% generative model
%==========================================================================
 
% first level
%--------------------------------------------------------------------------
G(1).f  = 'spm_fx_adem_write';
G(1).g  = 'spm_gx_adem_write';
G(1).x  = [pi/2; pi/2; 0; 0];                 % physical states
G(1).V  = exp(16);                            % error precision
G(1).W  = exp(16);                            % error precision
G(1).U  = sparse(1:4,1:4,exp(4),8,8);         % restriction
 
% second level
%--------------------------------------------------------------------------
G(2).v  = [0; 0];                             % exogenous forces
G(2).a  = [0; 0];                             % action forces
G(2).V  = exp(16);
 
 
% generate and invert
%==========================================================================
N       = 256;                                % length of data sequence
t       = (1:N)*8;
DEM.G   = G;
DEM.M   = M;
DEM.C   = sparse(2,N);
DEM     = spm_ADEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU,DEM.pU)
 
subplot(2,2,3)
spm_dem_reach_movie(DEM)
title('click on finger for movie','FontSize',16)
 
 
% Action-observation
%==========================================================================
ADEM          = DEM;
ADEM.C        = DEM.qU.a{2};
ADEM.M(1).V   = exp([-8 -8 -8 -8  8  8  8  8]);  % remove proprioception
ADEM.G(1).U   = 0;                               % remove motor control
ADEM.M(1).x.x = DEM.M(1).x.x;
ADEM.M(1).x.a = DEM.M(1).x.a*0;
ADEM          = spm_ADEM(ADEM);
 
spm_figure('GetWin','Figure 2');
spm_DEM_qU(ADEM.qU,ADEM.pU)
 
subplot(2,2,3)
spm_dem_reach_movie(ADEM)
title('click on finger for movie','FontSize',16)
 
 
% Hidden (intended) states lead actual states
%==========================================================================
spm_figure('GetWin','Figure 3');
 
n  = 64;
p  = exp(DEM.qU.x{1}(5:end,:));
qx = P(end,:)*(p*diag(1./sum(p)));
px = DEM.pU.v{1}(end,:);
 
subplot(2,1,1)
plot(t,px,t,qx,'-.')
legend({'attained','intended'})
title('intended and attained position','FontSize',16)
xlabel('time')
ylabel('vertical displacement')
box off
 
qx = qx - mean(qx);
px = px - mean(px);
 
subplot(2,2,3)
xc = xcorr(qx,px,n,'coef');
plot([0:(n + n)] - n,xc,[0 0],[-1 1],'k-.')
title('action','FontSize',16)
xlabel('time')
ylabel('cross-correlation')
axis square tight
 
 
% repeat for inferred hidden (intended) states
%--------------------------------------------------------------------------
p  = exp(ADEM.qU.x{1}(5:end,:));
qx = P(end,:)*(p*diag(1./sum(p)));
px = ADEM.pU.v{1}(end,:);
 
qx = qx - mean(qx);
px = px - mean(px);
 
subplot(2,2,4)
xca = xcorr(qx,px,n,'coef');
plot((0:(n + n)) - n,xca,(0:(n + n)) - n,xc,':',[0 0],[-1 1],'k-.')
title('observation','FontSize',16)
xlabel('time')
ylabel('cross-correlation')
axis square tight



% coherence analysis
%==========================================================================
spm_figure('GetWin','Figure 4');

T  = 16;
xe = DEM.qU.w{1}(5:end,T:end);      % prediction error on attractor states
ve = DEM.qU.z{1}(4:8,T:end);        % prediction error on visual states
qx = full(sum(xe));                 % synthetic LFP
qv = full(sum(ve));                 % synthetic LFP
px = DEM.pU.v{1}(1,T:end);          % peripheral (plant motion)

% coherence analysis
%--------------------------------------------------------------------------
qx = spm_detrend(qx',2);
px = spm_detrend(px',2);

[Cqp,Hz] = mscohere(qx,px,64,48,[],1/0.008);
[Cxv,Hz] = mscohere(qx,qv,64,48,[],1/0.008);


subplot(2,1,1)
plot(t(T:end),xe,'r',t(T:end),ve,'b');  hold on
plot(t(T:end),px/8,'g','LineWidth',4); hold off
title('central and peripheral responses','FontSize',16)
xlabel('time (ms)')
ylabel('activity (au)')
axis tight

subplot(2,2,3)
plot(Hz,Cqp)
title('central-peripheral coherence','FontSize',16)
xlabel('frequency (Hz)')
ylabel('coherence')
axis square
axis([0 32 0 1])

subplot(2,2,4)
plot(Hz,Cxv)
title('central-central coherence','FontSize',16)
xlabel('frequency (Hz)')
ylabel('coherence')
axis square
axis([0 32 0 1])

 
% functional correlates
%==========================================================================
spm_figure('GetWin','Figure 5');
 
% under action
%--------------------------------------------------------------------------
x    = DEM.pU.v{1}(7,:);
y    = DEM.pU.v{1}(8,:);
 
u    = DEM.qU.x{1}(8,:);
i    = find(u > 2);
 
subplot(2,2,1)
plot(x,y,'LineWidth',4,'Color',[1 1 1]*.8), hold on
plot(x(i),y(i),'r.','MarkerSize',24),       hold off
title('action','FontSize',16)
xlabel('position (x)')
ylabel('position (y)')
axis square
 
x  = x + ([1:N] - N)/N;
 
subplot(2,2,3)
plot(x,y,'LineWidth',4,'Color',[1 1 1]*.8), hold on
plot(x(i),y(i),'r.','MarkerSize',24),       hold off
title('action','FontSize',16)
xlabel('position (x)')
ylabel('position (y)')
axis square ij
 
 
% and observation
%--------------------------------------------------------------------------
x    = ADEM.pU.v{1}(7,:);
y    = ADEM.pU.v{1}(8,:);
 
u    = ADEM.qU.x{1}(8,:);
i    = find(u > 2);
 
subplot(2,2,2)
plot(x,y,'LineWidth',4,'Color',[1 1 1]*.8), hold on
plot(x(i),y(i),'r.','MarkerSize',24),       hold off
title('observation','FontSize',16)
xlabel('position (x)')
ylabel('position (y)')
axis square
 
x  = x + ([1:N] - N)/N;
 
subplot(2,2,4)
plot(x,y,'LineWidth',4,'Color',[1 1 1]*.8), hold on
plot(x(i),y(i),'r.','MarkerSize',24),       hold off
title('observation','FontSize',16)
xlabel('position (x)')
ylabel('position (y)')
axis square ij
 
 
% Correlations between action and observation responses
%==========================================================================
spm_figure('GetWin','Figure 6');
 
qx =  DEM.qU.x{1};
px = ADEM.qU.x{1};
 
subplot(2,1,1)
bar3(spm_en(px',0)'*spm_en(qx',0))
title('correlations','FontSize',16)
xlabel('hidden unit (observation)')
ylabel('hidden unit (action)')
grid off
 
% last hidden state
%--------------------------------------------------------------------------
qx =  DEM.qU.x{1}(end,:);
px = ADEM.qU.x{1}(end,:);
 
subplot(2,2,3)
plot(px,qx,'.')
title('correlations','FontSize',16)
xlabel('hidden unit (observation)')
ylabel('hidden unit (action)')
box off
axis square
 
qx = spm_en(qx',0);
px = spm_en(px',0);
 
subplot(2,2,4)
xc = xcorr(qx,px,n,'coef');
plot([0:(n + n)] - n,xc,[0 0],[-1 1],'k-.')
title('action','FontSize',16)
xlabel('time')
ylabel('cross-correlation')
axis square tight
 
 
 
% Simulated ERP
%==========================================================================
 
% simulate deviant
%--------------------------------------------------------------------------
T          = 132;
C          = DEM.qU.a{2};
C(:,T:end) = -C(:,T:end);
NDEM       = ADEM;
NDEM.C     = C;
NDEM       = spm_ADEM(NDEM);
 
 
% Graphics 
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7');
t   = ([1:N] - T)*8;
i   = find(t > -256 & t < 512);
 
subplot(3,2,1)
spm_dem_reach_movie(ADEM)
title('observation','FontSize',16)
 
subplot(3,2,2)
spm_dem_reach_movie(NDEM)
title('violation','FontSize',16)
 

% Proprioceptive prediction error 
%--------------------------------------------------------------------------
subplot(3,2,4)
plot(t(i),NDEM.qU.z{1}(1:4,i)'), hold on
plot(t(i),ADEM.qU.z{1}(1:4,i)',':'), hold off
title('proprioceptive error','FontSize',16)
xlabel('time')
ylabel('prediction error')
axis square
a = axis;
 
subplot(3,2,3)
plot(t(i),ADEM.qU.z{1}(1:4,i)')
title('proprioceptive error','FontSize',16)
xlabel('time')
ylabel('prediction error')
axis square
axis(a);
 
% Hidden state prediction error 
%--------------------------------------------------------------------------
subplot(3,2,6)
plot(t(i),NDEM.qU.w{1}(:,i)'), hold on
plot(t(i),ADEM.qU.w{1}(:,i)',':'), hold off
title('error on hidden states','FontSize',16)
xlabel('time')
ylabel('prediction error')
axis square
a = axis;
 
subplot(3,2,5)
plot(t(i),ADEM.qU.w{1}(:,i)')
title('error on hidden states','FontSize',16)
xlabel('time')
ylabel('prediction error')
axis square
axis(a);
