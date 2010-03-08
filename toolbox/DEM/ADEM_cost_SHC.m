% This demo illustrates the use of priors on the motion of hidden states as
% polices. It simulates exploration and exploitation using radial basis
% function attractors and auto-vitiative (self-destroying) attractors
% as the basis of the prior. These dynamics enforce exploration, under
% active inference. In turn, this exploration enables perceptual learning
% to associate attractors with changes in physiological states (cf,
% rewards). This can be exploited to by formal priors to ensure regions of
% physiological state-space are avoided.
% We look at this scheme using simulated pathology; first, we simulate a
% (neurodegenerative) reduction in log-precision (cf Parkinson’s disease) on
% the motion of physical states.  This results in active inference with
% a loss of precise volitional movement and subsequent failure to optimise 
% physiological states. Second, we look at the effects of precision on 
% learning by increasing log-precision (cf, Addition) on the motion of
% physiological states. This results in a failure to learn and, again, 
% subsequent failure to optimise physiological states.
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: ADEM_cost_SHC.m 3757 2010-03-08 11:41:53Z guillaume $
 
 
% switch for demo
%--------------------------------------------------------------------------
clear
DEMO = 1;
 
% location and radius of attractors (A)
%--------------------------------------------------------------------------
global A
A.x  = [1 -1 -1  1;
        2 -2  2 -2];
A.d  = 1/4;              % s.d. of Gaussian radial basis function
A.u  = 1/8;              % threshold to induce rapid decay of x.a
A.q  = [1 2];            % indices of locations that increase x.q
 
% generative process
%==========================================================================
 
% level 1
%--------------------------------------------------------------------------
G(1).x.x = [0;0];
G(1).x.v = [1;0];
G(1).x.q = [1;1/2];
 
G(1).f   = inline('spm_cost_SHC_fxa(x,v,a,P)','x','v','a','P');
G(1).g   = inline('[x.x; x.v; x.q]','x','v','a','P');
G(1).V   = exp(16);                          % error precision
G(1).W   = exp(16);                          % error precision
 
% level 2
%--------------------------------------------------------------------------
G(2).a   = sparse(2,1);                      % action
G(2).v   = 0;                                % inputs
G(2).V   = exp(16);
G        = spm_ADEM_M_set(G);
 
 
% generative model
%==========================================================================
 
% level 1
%--------------------------------------------------------------------------
M(1).x.x = G(1).x.x;
M(1).x.v = G(1).x.v;
M(1).x.q = G(1).x.q;
M(1).x.a = randn(4,1)/8;
 
M(1).f   = inline('spm_cost_SHC_fx(x,v,P)','x','v','P');
M(1).g   = inline('[x.x; x.v; x.q]','x','v','P');
M(1).pE  = speye(4,2);
 
M(1).V   = exp(8);
M(1).W   = diag(exp([[1 1 1 1]*4 [1 1]*6 [1 1 1 1]*4]));
 
% level 2 (no exogenous inputs in this simulation)
%--------------------------------------------------------------------------
M(2).v   = 0;
M(2).V   = exp(16);
 
% Integrate: active inference
%==========================================================================
M(1).E.nE  = 1;
M(1).E.n   = 4;
 
N      = 128;
DEM.U  = sparse(N,1);
DEM.C  = sparse(N,1);
DEM.G  = G;
DEM.M  = M;
 
if DEMO
    load DEM_addiction
else
    DEM = spm_ADEM(DEM);
end
 
 
% show behavior
%==========================================================================
 
% overview
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(DEM.qU)
 
subplot(2,2,3)
spm_cost_SHC_path(DEM.pU,A)
 
 
% a closer look at physiology
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2'); clf
 
subplot(2,1,1)
t    = 1:N;
plot(t,DEM.qU.x{1}(5:6,:)',t,t*0 + A.u,':'), hold on
plot(t,DEM.qU.x{1}(7:10,:)','-.'),           hold off
xlabel('time','Fontsize',14)
title('internal states','FontSize',16)
axis square, box off, set(gca,'XLim',[1 N])
 
subplot(2,2,3)
plot(t,DEM.qU.x{1}(5:6,:)',t,t*0 + A.u,':')
xlabel('time','Fontsize',14)
title('physiological states','FontSize',16)
axis square, box off, set(gca,'XLim',[1 N])
 
subplot(2,2,4)
plot(DEM.pU.x{1}(5,:),DEM.pU.x{1}(6,:)),               hold on
plot([A.u 1],[A.u A.u],'-.r',[A.u A.u],[A.u 1],'-.r'), hold off
xlabel('time','Fontsize',14)
title('physiological states','FontSize',16)
axis square, box off, axis([-.1 1.2 -.1 1.2])
 
 
 
% look at the effect of precision on inference (cf Parkinson’s disease)
%==========================================================================
 
% initialize states
%--------------------------------------------------------------------------
DEM = spm_ADEM_update(DEM);
 
% simulate exposure with different log-precisions on physical motion
%--------------------------------------------------------------------------
WP    = [4 2 0];
if ~DEMO
    for i = 1:length(WP)
        DEM_P{i}        = DEM;
        W               = [[1 1 1 1]*WP(i) [1 1]*6 [1 1 1 1]*4];
        DEM_P{i}.M(1).W = diag(exp(W));
        DEM_P{i}        = spm_ADEM(DEM_P{i});
    end
end
 
% Graphics - path and physiology
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3'); clf
for i = 1:3
    
    % path
    %----------------------------------------------------------------------
    spm_figure('GetWin','Figure 3');
    subplot(3,2,i*2 - 1)
    spm_cost_SHC_path(DEM_P{i}.pU,A)
    
    % and physiology
    %----------------------------------------------------------------------
    subplot(3,2,i*2)
    plot(t,DEM_P{i}.pU.x{1}(5:6,:)',t,t*0 + A.u,':'), hold on
    plot(t,sum(DEM_P{i}.pU.x{1}(1:2,:).^2)/8,'m:'),    hold off
    xlabel('time','Fontsize',12)
    title(sprintf('%s (%1.0f)','physiological states',WP(i)),'FontSize',16)
    axis square, box off, set(gca,'XLim',[1 N])
    drawnow
    
end
 
% Graphics - action and prediction error
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4'); clf
for i = 1:3
    
    % path
    %----------------------------------------------------------------------
    subplot(3,2,2*i - 1)
    plot(t,DEM_P{i}.qU.a{2})
    xlabel('time','FontSize',12)
    title('action','FontSize',16)
    axis square, box off, set(gca,'XLim',[1 N])
    if i == 1; a = axis; end, axis(a)
    
    % and prediction error
    %----------------------------------------------------------------------
    subplot(3,2,2*i)
    plot(t,DEM_P{i}.qU.z{1}(1:4,:))
    xlabel('time','FontSize',12)
    title('sensory error','FontSize',16)
    axis square, box off, set(gca,'XLim',[1 N])
    if i == 1; aa = axis; end, axis(aa)
    drawnow
end
 
% Learning (at different levels of precision on physiological motion)
%==========================================================================
 
% Switch attractor for second physiological state and log-precision levels
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM'); clf
A.q  = [1 3];
WA   = [4 8 12];
 
if ~DEMO
    
    for i = 1:length(WA)
 
        % Enable learning
        %------------------------------------------------------------------
        A.u              = 0;
        DEM_L{i}         = DEM;
        DEM_L{i}.M(1).pC = exp(8);
        W                = [[1 1 1 1]*4 [1 1]*WA(i) [1 1 1 1]*4];
        DEM_L{i}.M(1).W  = diag(exp(W));
        DEM_L{i}         = spm_ADEM(DEM_L{i});
 
 
        % replace prior with posterior and re-expose
        %------------------------------------------------------------------
        A.u              = 1/8;
        DEM_D{i}         = DEM_L{i};
        DEM_D{i}.M(1).pE = DEM_D{i}.qP.P{1};
        DEM_D{i}.M(1).pC = [];
        DEM_D{i}         = spm_ADEM(DEM_D{i});
 
    end
    
    % save DEM structures for future demonstrations
    %----------------------------------------------------------------------
    try
        save -v7 DEM_addiction DEM DEM_P DEM_L DEM_D
    catch
        save DEM_addiction DEM DEM_P DEM_L DEM_D
    end
    
end
 
 
% Graphics - optimal learning
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 5'); clf
spm_DEM_qP(DEM_L{1}.qP)
 
spm_figure('GetWin','Figure 5');
subplot(2,2,3)
spm_cost_SHC_path(DEM.pU,A)
title('Before','Fontsize',16)
subplot(2,2,4)
spm_cost_SHC_path(DEM_D{1}.pU,A)
title('After','Fontsize',16)
drawnow
 
 
% Comparison of learning over levels of log-precisions (cf. Addiction)
%==========================================================================
spm_figure('GetWin','Figure 6'); clf
 
for i = 1:3
    
    % path
    %----------------------------------------------------------------------
    subplot(3,2,i*2 - 1)
    spm_cost_SHC_path(DEM_D{i}.pU,A)
    title('After','FontSize',16)
    
    % and physiology
    %----------------------------------------------------------------------
    subplot(3,2,i*2)
    spm_plot_ci(DEM_L{i}.qP.P{1}(:),DEM_L{i}.qP.C)
    xlabel('parameter','FontSize',12)
    title(sprintf('%s (%1.0f)','learning',WA(i)),'FontSize',16)
    axis square, box off
    
end
drawnow

% and physiology
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 7'); clf
for i = 1:3
 
    % predicted motion of second physiological state
    %----------------------------------------------------------------------
    for j = 1:N
        x    = DEM_L{i}.qU.x{1}(:,j);
        f    = spm_cost_SHC_fx(spm_unvec(x,M(1).x),M(2).v,M(1).pE);
        p(j) = DEM_L{i}.qU.w{1}(6,j) + f.q(2);
        x    = DEM_L{i}.pU.x{1}(:,j);
        f    = spm_cost_SHC_fxa(spm_unvec(x,G(1).x),M(2).v,G(2).a,G(1).pE);
        q(j) = f.q(2);
    end
    subplot(3,2,2*i - 1)
    plot(t,p,t,q,'b:')
    xlabel('time','FontSize',12)
    title(sprintf('%s (%1.0f)','predicted motion',WA(i)),'FontSize',16)
    axis square, box off, set(gca,'XLim',[1 N])
    if i == 1; a = axis; end, axis(a)
 
    subplot(3,2,2*i)
    plot(t,DEM_L{i}.qU.w{1}(6,:))
    xlabel('time','FontSize',12)
    title('prediction error','FontSize',16)
    axis square, box off, set(gca,'XLim',[1 N])
    if i == 1; aa = axis; end, axis(aa)
 
end
drawnow