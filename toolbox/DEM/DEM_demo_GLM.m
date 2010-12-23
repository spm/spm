% Demo comparing DEM and ReML (restricted maximum likelihood) under a simple
% general linear model (GLM).  Slight differences in the hyperpriors of both
% schemes make this an interesting exercise.  Note that ReML uses a
% covariance hyper-parameterisation; whereas DEM uses precision
% hyperparameters.  This demo uses a non-hierarchical GLM and switches the
% roles of parameters and causes to illustrate their equivalence under a DEM
% inversion.
 
% Classical (OLS) estimation of states');
%==========================================================================
 
% specify parameters
%--------------------------------------------------------------------------
N     = 16;                                        % length of data sequence
X     = randn(8,2);
 
% generate data
%--------------------------------------------------------------------------
DEM   = spm_DEM_generate(spm_DEM_M('GLM',X),N,{},{1});
 
% DEM estimation
%==========================================================================
DEM   = spm_DEM(DEM);
qU    = DEM.qU;
qP    = DEM.qP;
qH    = DEM.qH;
 
% Classical (OLS) estimation of states
%==========================================================================
[Ce,h,W,L] = spm_reml(DEM.Y*DEM.Y',X,DEM.M(1).Q,N);
V          = pinv(X)*DEM.Y;
C          = h*pinv(X)*DEM.M(1).Q{1}*pinv(X)';
 
% compare estimators - hyperparameters
%--------------------------------------------------------------------------
subplot(3,2,1)
c = sqrt(qH.C)*spm_invNcdf(1 - 0.05);
bar(full(exp(-qH.h{1})),'c')
line([1 1],exp(-([-1 1]*c + qH.h{1})),'LineWidth',4,'Color','r');
axis square; title('DEM hyperparameter estimates')
a = axis;

subplot(3,2,2)
c = sqrt(inv(W))*spm_invNcdf(1 - 0.05);
bar(full(h),'c')
line([1 1],exp(-([-1 1]*c - log(h))),'LineWidth',4,'Color','r');
axis square; title('ReML hyperparameter estimates')
axis(a)
 
% and causes
%--------------------------------------------------------------------------
subplot(3,2,3)
bar(qU.v{2})
axis square; title('DEM parameter estimates')
a = axis;

subplot(3,2,4)
bar(V)
axis square; title('ReML parameter estimates')
axis(a)
 
subplot(3,2,5)
imagesc(qU.C{1})
axis square; title('DEM Cq')

subplot(3,2,6)
imagesc(C)
axis square; title('ReML Cq')

 
% transpose model: states->parameters to test E-Step
%========================================================================== 
q = questdlg('transpose (states->parameters) GLM to test E-Step');
 
if ~strcmp(q,'Yes'), return, end
%==========================================================================
M         = spm_DEM_M('GLM',sparse(N,size(X,2)));
M(2).V    = speye(2,2)*exp(8);
M(1).pC   = speye(size(X,2)*N)*exp(16);
M(1).E.nM = 1;

Y       = DEM.Y';
TEM.M   = M;
TEM.Y   = Y;
TEM.U   = X';
 
% DEM estimation hyperparameters
%==========================================================================
TEM     = spm_DEM(TEM);
qU      = TEM.qU;
qP      = TEM.qP;
qH      = TEM.qH;
 
% compare estimators - hyperparameters
%--------------------------------------------------------------------------
subplot(3,2,1)
c = sqrt(qH.C)*spm_invNcdf(1 - 0.05);
bar(full(exp(-qH.h{1})),'c')
line([1 1],exp(-([-1 1]*c + qH.h{1})),'LineWidth',4,'Color','r');
axis square; title('DEM hyperparameter estimates')
a = axis;

subplot(3,2,2)
c = sqrt(inv(W))*spm_invNcdf(1 - 0.05);
bar(full(h),'c')
line([1 1],exp(-([-1 1]*c - log(h))),'LineWidth',4,'Color','r');
axis square; title('ReML hyperparameter estimates')
axis(a)

% and causes (i.e., parameters)
%--------------------------------------------------------------------------
subplot(3,2,3)
bar(qP.P{1}')
axis square; title('DEM parameter estimates')
a = axis;

subplot(3,2,4)
bar(V)
axis square; title('ReML parameter estimates')
axis(a)
 
subplot(3,2,5)
imagesc(qP.C)
axis square; title('DEM Cq')

subplot(3,2,6)
imagesc(kron(C,eye(N)))
axis square; title('ReML Cq')
