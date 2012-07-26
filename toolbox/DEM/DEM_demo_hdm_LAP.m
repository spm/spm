function DEM_demo_hdm_LAP
% Demo for Hemodynamic deconvolution: Cross-validation of Laplace scheme
%__________________________________________________________________________
% This demonstration compares generalised filtering and DEM in the context
% of a nonlinear convolution model using empirical data. These are the data
% used to illustrate hemodynamic deconvolution. We have deliberately made
% the problem difficult here to highlight the ability of Generalised
% filtering to accumulate evidence to optimise in parameters and hyper-
% parameters, which allows it to outperform DEM (although it does not
% find visual motion effects with 90% confidence)
%__________________________________________________________________________
% Copyright (C) 2010 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: DEM_demo_hdm_LAP.m 4804 2012-07-26 13:14:18Z karl $

% set-up
%--------------------------------------------------------------------------
global dt

spm_figure('GetWin','DEM');
load HDM

% generative [likelihood] model 'HDM'
%==========================================================================
dt      = Y.dt;
T       = 1:256;

% level 1
%--------------------------------------------------------------------------
[pE pC] = spm_hdm_priors(3,2);
M(1).x  = [0 0 0 0]';
M(1).g  = 'spm_gx_hdm';
M(1).f  = 'spm_fx_hdm';
M(1).pE = pE;
M(1).pC = pC;
M(1).W  = exp(8);
M(1).V  = exp(4);

% level 2
%--------------------------------------------------------------------------
M(2).v  = [0 0 0]';
M(2).V  = exp(2);

M(1).E.linear = 1;
M(1).E.n  = 4;
M(1).E.nD = 1;
M(1).E.nE = 16;
M(1).E.nN = 16;
M(1).E.s  = 1/2;

% Decimate U.u from micro-time
% -------------------------------------------------------------------------
Dy  = spm_dctmtx(size(Y.y,1),size(Y.y,1));
Du  = spm_dctmtx(size(U.u,1),size(Y.y,1));
Dy  = Dy*sqrt(size(Y.y,1)/size(U.u,1));
u   = Dy*(Du'*U.u);
U   = spm_detrend(u(T,:))';


% Emprical data
%--------------------------------------------------------------------------
DEM.M = M;
DEM.U = U;
DEM.Y = Y.y(T,:)'/4;

% DEM estimation
%==========================================================================
DEM        = spm_DEM(DEM);
LAP        = spm_LAP(DEM);



% report states and parameter estimates
%==========================================================================

% LAP
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 1');
spm_DEM_qU(LAP.qU)

% Coupling parameters of interest
%--------------------------------------------------------------------------
subplot(2,2,4)
qP    = LAP.qP.P{1}(7:end);
bar(qP,'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8)
cq    = 1.64*sqrt(diag(LAP.qP.C(7:end,7:end)));
hold on
for i = 1:length(qP)
    plot([i i], qP(i) + [-1 1]*cq(i),'LineWidth',8,'color','r')
end
hold off
axis square
set(gca,'XTickLabel',{'vision','motion','attention'})
title({'parameters';'Generalised filtering (Laplace)'},'Fontsize',16)
a   = axis;

% DEM
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 2');
spm_DEM_qU(DEM.qU)


% Coupling parameters of interest
%--------------------------------------------------------------------------
subplot(2,2,4)
qP    = DEM.qP.P{1}(7:end);
bar(qP,'Edgecolor',[1 1 1]/2,'Facecolor',[1 1 1]*.8)
cq    = 1.64*sqrt(diag(DEM.qP.C(7:end,7:end)));
hold on
for i = 1:length(qP)
    plot([i i], qP(i) + [-1 1]*cq(i),'LineWidth',8,'color','r')
end
hold off
axis square
set(gca,'XTickLabel',{'vision','motion','attention'})
title({'parameters';'Mean-field approximation (DEM)'},'Fontsize',16)
axis(a)


% Log-evidence
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 3');

subplot(2,1,1)
nL   = length(LAP.F);
nD   = length(DEM.F);
plot(1:nL,LAP.F,1:nD,DEM.F)
axis square
legend('LAP (F)','DEM(F)')
title('log-evidence','FontSize',16)


return


% and a more detailed look
%--------------------------------------------------------------------------
spm_figure('GetWin','Figure 4: Hemodynamics');


t   = 1:128;
subplot(2,1,1)
hold on
bar(full(LAP.U(2,t)*8),'FaceColor',[1 1 1]*.8,'EdgeColor',[1 1 1]*.8)
plot(t,exp(LAP.qU.x{1}(:,t)))
set(gca,'YLim',[-0.1 1.6])
xlabel('time (bins)','Fontsize',12)
title('hidden states','Fontsize',16)
legend({'visual stimulation','signal','flow','volume','dHb'})
hold off

% (mixture of) causes
%--------------------------------------------------------------------------
qP  = LAP.qP.P{1}(7:end);

subplot(2,1,2)
hold on
plot(t,qP'*LAP.qU.v{2}(:,t))
a = axis;
bar(full(LAP.U(2,t)),'FaceColor',[1 1 1]*.8,'EdgeColor',[1 1 1]*.8)
plot(t,qP'*LAP.qU.v{2}(:,t))
axis(a)
xlabel('time (bins)','Fontsize',12)
title('neuronal causes','Fontsize',16)
hold off


return


% Notes for simulating responses (to examine dependency on
% hemodynamic parameters quantitatively)
%==========================================================================

% get causes (i.e. experimental inputs)
%--------------------------------------------------------------------------
T  = 128;
P  = [
    0.9874
    0.3501
    1.6845
    0.3452
    0.3564
    -0.1719
    0.0000
    0.2
    0.0000];

P  = [
    0.6500
    0.4100
    2.0000
    0.3200
    0.3400
    -.0
    0
    .16
    0];


% Simulate data
%--------------------------------------------------------------------------
DEM    = spm_DEM_generate(M,U(:,1:T),{P},{4,8},{6});
spm_DEM_qU(DEM.pU)

