% demo for Hemodynamic deconvolution: Validation of Laplace scheme
%__________________________________________________________________________

% set-up
%==========================================================================
clear
spm_figure('GetWin','DEM');
load HDM
global dt


% generative [likelihood] model 'HDM'
%==========================================================================
dt      = Y.dt;
T       = 1:256;

% level 1
%--------------------------------------------------------------------------
[pE pC] = spm_hdm_priors(3,1);
pC(6,6) = 0;
M(1).x  = [0 0 0 0]';
M(1).g  = 'spm_gx_hdm';
M(1).f  = 'spm_fx_hdm';
M(1).pE = pE;
M(1).pC = pC;
M(1).hE = 2;
M(1).hC = 1;
M(1).gE = 2;
M(1).gC = 1;

% level 2
%--------------------------------------------------------------------------
M(2).v  = [0 0 0]';
M(2).V  = 1;

M(1).E.linear = 1;
M(1).E.n  = 4;
M(1).E.nD = 1;
M(1).E.nN = 16;
M(1).E.s  = 1/4;

% get causes (i.e. experimental inputs)
%--------------------------------------------------------------------------
t  = fix(linspace(1,length(U.u),360));
U  = U.u(t(T),:)';

% Emprical data
%--------------------------------------------------------------------------
DEM.M  = M;
DEM.U  = U;
DEM.Y  = 1 + Y.y(T,:)'/4;

% DEM estimation
%==========================================================================
LAP    = spm_LAP(DEM);
DEM    = spm_DEM(DEM);



% report states and parameter esimates
%==========================================================================

% DEM
%--------------------------------------------------------------------------
spm_figure('GetWin','DEM');
spm_DEM_qU(DEM.qU)

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
title('parameters','Fontsize',16)


% DEM
%--------------------------------------------------------------------------
spm_figure('GetWin','Laplace');
spm_DEM_qU(LAP.qU)

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
title('parameters','Fontsize',16)


% Log-evidence
%--------------------------------------------------------------------------
spm_figure('GetWin','SI');
clf

subplot(2,1,1)
nL   = length(LAP.F);
nD   = length(DEM.F);
plot(1:nL,LAP.S,1:nD,DEM.S)
axis square
legend('LAP (S)','DEM(S)')
title('log-evidence ','FontSize',16)


return


% and a more detailed look
%--------------------------------------------------------------------------
t = 1:128;
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



% Simulated respsones
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


