function ALAP_demo_attenuation
% This demonstration illustrates context or state-dependent precision (i.e.
% attention), which is necessary to disambiguate between sensations
% caused exogenously and self-generated sensations. In brief, it is
% necessary to attend away from the sensory consequences of action to
% preclude sensory evidence overriding the prior beliefs that cause
% movement. This necessarily reduced the confidence in self-generated
% sensations and provides a simple (Bayes-optimal) explanation for sensory
% attenuation - in terms of the attention of sensory precision. We
% illustrate this in the setting of the force matching illusion and go on
% to show that increasing the conviction in (precision of) prior beliefs
% abolishes sensory attenuation at the expense of false (delusional) 
% posterior beliefs about antagonistic external forces.
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ALAP_demo_attenuation.m 4826 2012-08-03 16:45:09Z karl $
 
% process (G) and model (M)
%==========================================================================
 
% set dimensions for generalised coordinates
%--------------------------------------------------------------------------
G(1).E.d        = 2;                   % approximation order
G(1).E.n        = 4;                   % embedding order
 
M(1).E.d        = 2;                   % approximation order
M(1).E.n        = 4;                   % embedding order
M(1).E.s        = 1/2;                 % embedding order
 
M(1).E.method.x = 1;                   % state-dependent noise
M(1).E.method.v = 1;                   % state-dependent noise
M(1).E.method.h = 0;                   % suppress optimisation
M(1).E.method.g = 0;                   % suppress optimisation
 
G(1).f  = inline('tanh(a) - x/4','x','v','a','P');
G(1).g  = inline('[x; v + x]','x','v','a','P');
G(1).x  = 0;                           % hidden state
G(1).v  = [0; 0];                      % hidden cause (sensory data)
G(1).V  = exp(8);                      % precision (noise)
G(1).W  = exp(8);                      % precision (states)
G(1).U  = [exp(0) 0];                  % precision (action)
 
 
% level 2; causes
%--------------------------------------------------------------------------
G(2).v  = 0;                           % exogenous  cause
G(2).a  = 0;                           % endogenous cause (action)
G(2).V  = exp(16);
 
 
% state-dependent precision (attentional bias) in generative model (M):
%--------------------------------------------------------------------------
M(1).f  = inline('v - x/4','x','v','P');
M(1).g  = inline('[x(1); sum(x)]','x','v','P');
M(1).x  = [0; 0];                      % hidden states
M(1).v  = [0; 0];                      % hidden causes
M(1).W  = exp(4);                      % precision (states)
M(1).ph = inline('[1; 1]*(8 - h*tanh(v(1) + x(1)))','x','v','h','M');
M(1).hE = 6;
 
 
% level 2; causes
%--------------------------------------------------------------------------
M(2).v  = [0; 0];                      % hidden cause
M(2).V  = [exp(6); exp(0)];
 
 
 
 
% Demonstration of the need for sensory attenuation
%==========================================================================
 
% hidden cause and prior expectations
%--------------------------------------------------------------------------
N      = 32;
C      = zeros(1,N);
U(1,:) = exp(-((1:N) - N/2).^2/(4.^2))*1;
U(2,:) = zeros(1,N);
 
% assemble model structure
%--------------------------------------------------------------------------
DEM.M = M;
DEM.G = G;
DEM.C = C;
DEM.U = U;
 
hE    = (-4:1:6);
for i = 1:length(hE)
    
    rng('default')
    
    LAP         = DEM;
    LAP.M(1).hE = hE(i);
    LAP         = spm_ALAP(LAP);
    
    % true and perceived force exerted (endogenously)
    %----------------------------------------------------------------------
    Px(i) = max(LAP.pU.x{1}(1,:));
    Qx(i) = max(LAP.qU.x{1}(1,:));
    
    
    % plot self-generated movement
    %----------------------------------------------------------------------
    if hE(i) == 0
        
        spm_figure('GetWin','Figure 1: Low attenuation');
        spm_DEM_qU(LAP.qU,LAP.pU)
        
    elseif hE(i) == 6
        
        spm_figure('GetWin','Figure 2: High attenuation');
        spm_DEM_qU(LAP.qU,LAP.pU)
        
    end
end
 
% adjust axes
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2: High attenuation');
subplot(2,2,2); spm_axis tight, a = axis;
subplot(2,2,1); axis(a);
spm_figure('GetWin','Figure 1: Low attenuation');
subplot(2,2,2); axis(a);
subplot(2,2,1); axis(a);
 
spm_figure('GetWin','Figure 2: High attenuation');
subplot(2,2,3); spm_axis tight, a = axis;
subplot(2,2,4); axis(a);
spm_figure('GetWin','Figure 1: Low attenuation');
subplot(2,2,3); axis(a);
subplot(2,2,4); axis(a);
 
 
 
% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');
 
subplot(2,1,1)
plot(hE,[Px; Qx])
axis square
xlabel('attenuation of sensory precicion','FontSize',12)
ylabel('true and perceived force exerted','FontSize',12)
legend({'true','perceived'})
title('Sensory attenuation and action','FontSize',16)
 
 
 
% Demonstration of sensory attenuation
%==========================================================================
 
% replay internal force as external force
%--------------------------------------------------------------------------
rng('default')
DEM.C  = [C LAP.pU.x{1}];
DEM.U  = [U sparse(2,N)];
DEM    = spm_ALAP(DEM);
 
% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4');
spm_DEM_qU(DEM.qU,DEM.pU)
 
 
 
% Force matching
%==========================================================================
 
% replay perceived (at 90% confidence) internal force
%--------------------------------------------------------------------------
for i = 1:N
    CI(i) = 1.694*sqrt(LAP.qU.S{i}(1,1));
end
 
DEM.C  = [C (LAP.qU.x{1}(1,:) - CI)];
DEM.U  = [U sparse(2,N)];
DEM    = spm_ALAP(DEM);
 
% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5');
spm_DEM_qU(DEM.qU,DEM.pU)
subplot(2,2,2); spm_axis tight, a = axis;
subplot(2,2,1); axis(a);
subplot(2,2,3); spm_axis tight, a = axis;
subplot(2,2,4); axis(a);
 
 
 
% Use a range of self-generated forces - with D = 0
%==========================================================================
F     = (1:4)/2;
for i = 1:length(F)
    
    DEM.C = C;
    DEM.U = U*F(i);
    DEM   = spm_ALAP(DEM);
    
    % self-generated and matched (inferred) force
    %----------------------------------------------------------------------
    [x j] = max(DEM.pU.x{1}(1,:));
    Sx(i) = x;
    Tx(i) = DEM.qU.x{1}(1,j) - 1.694*sqrt(DEM.qU.S{j}(1,1));
    
end
 
% repeat with (delusional) precision D = 2
%--------------------------------------------------------------------------
D           = 2;
DEM.M(1).hE = 6 - D;
DEM.M(1).W  = exp(4 + D);
DEM.M(2).V  = [exp(6 + D); exp(0)];
 
 
for i = 1:length(F)
    
    DEM.C = C;
    DEM.U = U*F(i);
    DEM   = spm_ALAP(DEM);
    
    % self-generated and matched (inferred) force
    %----------------------------------------------------------------------
    [x j] = max(DEM.pU.x{1}(1,:));
    sx(i) = x;
    tx(i) = DEM.qU.x{1}(1,j) - 1.694*sqrt(DEM.qU.S{j}(1,1));
    
end
 
% plot results of force matching
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 6');
 
subplot(2,1,1)
plot(Tx,Sx,'-'); hold on,
plot(Tx,Sx,'.','MarkerSize',32), hold on
plot(tx,sx,'-.'); hold on,
plot(tx,sx,'.','MarkerSize',16), hold on
plot([0 3],[0 3],':'); hold off
 
axis square
xlabel('external (perceptually matched) force','FontSize',12)
ylabel('self-generated force','FontSize',12)
title('Force matching illusion','FontSize',16)
 
% Illustrate false inference with delusional precision
%==========================================================================
D           = 4;
DEM.M(1).hE = 6 - D;
DEM.M(1).W  = exp(4 + D);
DEM.M(2).V  = [exp(6 + D); exp(0)];
 
 
% replay perceived (at 90% confidence) internal force
%--------------------------------------------------------------------------
DEM.C = C;
DEM.U = U*2;
DEM   = spm_ALAP(DEM);
 
for i = 1:N
    CI(i) = 1.694*sqrt(DEM.qU.S{i}(1,1));
end
 
DEM.C = [C (DEM.qU.x{1}(1,:) - CI)];
DEM.U = [U*2 sparse(2,N)];
DEM   = spm_ALAP(DEM);
 
 
% plot false inference under high D
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7: False (delusional) inference');
spm_DEM_qU(DEM.qU,DEM.pU)
