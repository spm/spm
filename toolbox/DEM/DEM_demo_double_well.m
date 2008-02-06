% DEMO comparing DEM with particle filtering in the context of a bimodal
% conditional density.  This demonstrates a shortcoming of DEM in that it
% fails to represent the true density.

 
% get nonlinear state-space model
%==========================================================================
Fgraph = spm_figure('GetWin','Graphics');
M      = spm_DEM_M('ssm');
 
% generate data (output)
%--------------------------------------------------------------------------
T      = 64;
U      = 8*sin(pi*[1:T]/16);
DEM    = spm_DEM_generate(M,U);
spm_DEM_qU(DEM.pU);
 
% EKF
%--------------------------------------------------------------------------
[kf_x] = spm_ekf(M,DEM.Y);
 
% PF
%--------------------------------------------------------------------------
[pf_x,P,Q,xQ] = spm_pf(M,DEM.Y);
 
% DEM
%--------------------------------------------------------------------------
DEM    = spm_DEM(DEM);
de_x   = DEM.qU.x{1};
tr_x   = DEM.pU.x{1};
 
 
% Graphical comparison
%--------------------------------------------------------------------------
figure(Fgraph)
t      = 1:T;
subplot(2,2,1)
plot(t,pf_x,t,kf_x,':',t,de_x,t,tr_x)
legend({'PF','EKF','DEM','true'})
title('hidden state')
axis([1 T -32 32])
axis square
 
subplot(2,2,2)
plot(t,tr_x - pf_x,t,tr_x - kf_x,':',t,tr_x - de_x)
legend({'PF','EKF','DEM'})
title('error')
axis([1 T -32 32])
axis square
 
% Sample density
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc([1:T],xQ,Q)
axis xy square
title('sample density')
