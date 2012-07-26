function DFP_demo_hdm
% demo for Hemodynamic deconvolution usinf variational filtering
%__________________________________________________________________________

% conventional analysis % [Ep,Cp,Ce,K0,K1,K2,M0,M1] = spm_nlsi(M,U,Y);
%==========================================================================
load HDM
global dt


% generative [likelihood] model 'HDM'
%==========================================================================
G(1).E.linear = 1;
G(1).E.s      = 1/2;
dt            = Y.dt;
T             = 1:128;


% level 1
%--------------------------------------------------------------------------
pE = [
    1.20
    0.31
    2.14
    0.36
    0.36
    0.84
    1];

G(1).x  = [0 0 0 0]';
G(1).g  = 'spm_gx_hdm';
G(1).f  = 'spm_fx_hdm';
G(1).pE = pE;
G(1).V  = speye(1,1)*exp(2);
G(1).W  = speye(4,4)*exp(8);

% level 2
%--------------------------------------------------------------------------
G(2).v  = 0;
G(2).V  = 1;

% get causes (i.e. experimental inputs)
%--------------------------------------------------------------------------
t      = fix(linspace(1,length(U.u),360));
DEM.M  = G;
DEM.Y  = Y.y(T,:)'/2;

% DEM estimation
%==========================================================================
DEM.M(1).E.N  = 16;
DEM.M(1).E.xp = 128;
DEM.M(1).E.xp = 128;

DFP    = spm_DFP(DEM);
DEM    = spm_DEM(DEM);

% states and parameter esimates
%==========================================================================
spm_DEM_qU(DEM.qU)


return

% and a more detailed look
%--------------------------------------------------------------------------
t      = fix(linspace(1,length(U.u),360));
U      = U.u(t(T),:)';

spm_figure('GetWin','Figure 2');
clf
subplot(2,1,1)
hold on
bar(full(U(2,:)*8),'FaceColor',[1 1 1]*.8,'EdgeColor',[1 1 1]*.9)
plot(T,exp(DFP.qU.x{1}))
set(gca,'YLim',[-0.1 1.6])
xlabel('time (bins)')
title('hidden states')
legend({'visual stimulation','signal','flow','volume','dHb'})
hold off

% causes
%--------------------------------------------------------------------------
subplot(2,1,2)
hold on
plot(T,DEM.qU.v{2})
a = axis;
bar(full(U(2,:)),'FaceColor',[1 1 1]*.8,'EdgeColor',[1 1 1]*.9)
plot(T,DEM.qU.v{2},'Linewidth',2)
axis(a)
xlabel('time (bins)')
title('neuronal causes')
hold off


