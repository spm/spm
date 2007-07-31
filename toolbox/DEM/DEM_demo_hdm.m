% demo for Hemodynamic deconvolution
%__________________________________________________________________________

% conventional analysis % [Ep,Cp,Ce,K0,K1,K2,M0,M1] = spm_nlsi(M,U,Y);
%==========================================================================
clear
Fgraph  = spm_figure('GetWin','Graphics');
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
[pE pC] = spm_hdm_priors(3,1);
G(1).x  = [0 0 0 0]';
G(1).g  = 'spm_gx_hdm';
G(1).f  = 'spm_fx_hdm';
G(1).pE = pE;
G(1).pC = pC;
G(1).Q  = speye(1,1);
G(1).R  = speye(4,4);
G(1).gE = 4;
G(1).gC = 1/32;

% level 2
%--------------------------------------------------------------------------
G(2).v  = [0 0 0]';
G(2).V  = speye(3,3);


% get causes (i.e. experimental inputs)
%--------------------------------------------------------------------------
t      = fix(linspace(1,length(U.u),360));
DEM.U  = U.u(t(T),:)';
DEM.M  = G;
DEM.Y  = Y.y(T,:)'/2;
DEM.X  = Y.X0(T,1)';

% DEM estimation
%==========================================================================
DEM    = spm_DEM(DEM);

% states and parameter esimates
%==========================================================================
spm_DEM_qU(DEM.qU)


subplot(2,2,4)
qP    = DEM.qP.P{1}(7:end);
bar(qP)
cq    = 1.64*sqrt(diag(DEM.qP.C(7:end,7:end)));
hold on
for i = 1:length(qP)
    plot([i i], qP(i) + [-1 1]*cq(i),'LineWidth',8,'color','r')
end
hold off
axis square
set(gca,'XTickLabel',{'vis','mot','att'})
title('parameters')


% and a more detailed look
%--------------------------------------------------------------------------
figure(Fgraph)
clf
t = 1:128;
subplot(2,1,1)
hold on
bar(full(DEM.U(2,t)*8),'FaceColor',[1 1 1]*.8,'EdgeColor',[1 1 1]*.9)
plot(t,exp(DEM.qU.x{1}(:,t)))
set(gca,'YLim',[-0.1 1.6])
xlabel('time (bins)')
title('hidden states')
legend({'visual stimulation','signal','flow','volume','dHb'})
hold off

% (mixture of) causes
%--------------------------------------------------------------------------
qP  = DEM.qP.P{1}(7:end);

subplot(2,1,2)
hold on
plot(t,qP'*DEM.qU.v{2}(:,t))
a = axis;
bar(full(DEM.U(2,t)),'FaceColor',[1 1 1]*.8,'EdgeColor',[1 1 1]*.9)
plot(t,qP'*DEM.qU.v{2}(:,t),'Linewidth',2)
axis(a)
xlabel('time (bins)')
title('neuronal causes')
hold off


