% DEM demo for a hierarchical linear model (MFX).  This inversion is
% cross-validated with restricted maximum likelihood using and parametric 
% empirical Bayes (spm_PEB). It uses a simple two level model that embodies
% empirical shrinkage priors on the first-level parameters (c.f.,
% DEM_demo_GLM, with no empirical priors)

% MFX design
%==========================================================================
X1    = kron(eye(4),kron(ones(2,1),eye(2)));
X2    = kron(ones(4,1),eye(2));
 
% specify parameters
%--------------------------------------------------------------------------
N     = 1;                                      % length of data sequence
h     = {log(16) log(8)};                       % precisions
 
% generate data
%--------------------------------------------------------------------------
DEM   = spm_DEM_generate(spm_DEM_M('HLM',X1,X2),N,{},h);
 
% DEM estimation
%==========================================================================
DEM.M(1).E.nE = 16;

[DEM] = spm_DEM(DEM);
qU    = DEM.qU;
qP    = DEM.qP;
qH    = DEM.qH;
 
% compare real and estimated factor and causes (using ReML and PEB)
%==========================================================================
Q1           = DEM.M(1).Q{1};
Q2           = DEM.M(2).Q{1};
[Ce,hr,W,Fr] = spm_reml_sc(DEM.Y*DEM.Y',X1*X2,{Q1; X1*Q2*X1'},N);
P{1}.X       = X1;
P{1}.C       = {Q1};
P{2}.X       = X2;
P{2}.C       = {Q2};
[C,P,Fp]     = spm_PEB(DEM.Y,P);
 
% graphics for variance component estimators
%--------------------------------------------------------------------------
subplot(3,1,1)
qh   = spm_vec(qH.h);
bar([exp(-[h{1}; h{2}]) [exp(-qh)] [C{1}.h; C{2}.h]])
legend({'true' 'DEM' 'ReML'})
title('ReML estimators')
c     = sqrt(diag(qH.C))*1.64;
for i = 1:length(h)
    line([i i],exp(-([-c(i) c(i)] + qh(i))),'color',[1 0 0],'LineWidth',2)
end
axis square
grid on
 
% check estimators
%--------------------------------------------------------------------------
subplot(3,2,3)
bar(qU.v{2})
axis square; title('DEM level-1')
a = axis;
subplot(3,2,4)
bar(C{2}.E)
axis square; title('VB level-1')
axis(a)
 
subplot(3,2,5)
bar(qU.v{3})
axis square; title('DEM level-2')
a = axis;
subplot(3,2,6)
bar(C{3}.E)
axis square; title('VB level-2')
axis(a)

