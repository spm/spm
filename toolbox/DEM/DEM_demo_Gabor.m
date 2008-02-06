% State-space demo routine simulating position invariance representations
% in the visual system.  The generative model predicts a one-dimensional
% Gabor patch that moves in the (one-dimensional) visual field. The
% inversion of this dynamic model can be viewed as deconvolving spatial and
% category attributed from a moving stimulus (or selective re-sampling of
% the input to recover the stimulus that can be represented.
%==========================================================================
clear M
Fgraph    = spm_figure('GetWin','Graphics');
 
% temporal correlations
%--------------------------------------------------------------------------
M(1).E.s  = 1/2;
 
 
 
% model specification - 1st level
%--------------------------------------------------------------------------
M(1).x  = [0;0;0];
M(1).f  = 'spm_fx_Gabor';
M(1).g  = 'spm_gx_Gabor';
M(1).V  = exp(8);
M(1).W  = exp(8);
 
% model specification - 2nd level
%--------------------------------------------------------------------------
M(2).v  = [0;0;0];
M(2).V  = 2;
 
% generate data (output)
%--------------------------------------------------------------------------
T       = 64;
U       = sparse(3,T);
U(1,:)  = 4*sin(pi*[1:T]/16);
U(3,:)  = 2*sin(pi*[1:T]/8);
DEM     = spm_DEM_generate(M,U);
spm_DEM_qU(DEM.pU)
 
subplot(2,2,4)
imagesc(DEM.Y)
axis square
xlabel time
 
% DEM
%==========================================================================
DEM        = spm_DEM(DEM);
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
 
% get prediction without displacement
%--------------------------------------------------------------------------
qx      = DEM.qU.x{1};
px      = DEM.pU.x{1};
qx(1,:) = qx(1,:) - px(1,:);
px(1,:) = px(1,:) - px(1,:);
for i = 1:T
    P(:,i) = spm_gx_Gabor(px(:,i),[],[]);
    Q(:,i) = spm_gx_Gabor(qx(:,i),[],[]);
end
 
% plot
%--------------------------------------------------------------------------
Fgraph    = spm_figure('GetWin','Graphics');
 
subplot(2,2,1)
imagesc(DEM.Y)
axis square
xlabel time
 
subplot(2,2,2)
imagesc(DEM.qU.v{1})
axis square
xlabel time
 
subplot(2,2,3)
imagesc(P)
axis square
xlabel time
 
subplot(2,2,4)
imagesc(Q)
axis square
xlabel time
