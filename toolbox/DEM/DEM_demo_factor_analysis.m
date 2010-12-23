% Demo for Probabilistic Factor Analysis; This uses a hierarchical model
% under the constraint that the causes have a deterministic and stochastic
% components.  The aim is to recover the true subspace of the real causes.
 
% non-hierarchical linear generative model (static)
%==========================================================================
n     = 8;
m     = 2;
M     = spm_DEM_M('FA',[n m]);
 
% create data
%==========================================================================
N     = 8;                                        % length of data sequence
X     = randn(size(M(1).pE));
DEM   = spm_DEM_generate(M,N,{X},{4});
 
% Initialise parameters
%--------------------------------------------------------------------------
DEM.class = 'FA';
DEM   = spm_dem_initialise(DEM);
 
% DEM estimation
%==========================================================================
DEM.M(1).E.nE = 16;
DEM   = spm_DEM(DEM);
 
% compare real and estimated factor and causes
%==========================================================================
 
% plot
%--------------------------------------------------------------------------
subplot(2,1,1)
v     = DEM.qU.v{2};
u     = DEM.pU.v{2};
plot(v'*pinv(full(v'))*u')
hold on
plot(u',':')
title('real and rotated causes')
axis square
grid on
