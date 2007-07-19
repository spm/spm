% DEMO comparing DEM with particle filtering in the context of a bimodal
% conditional density.  This demonstrates a shortcoming of DEM in that it
% fails to represent the true density
%==========================================================================
clf
clear M
 
% phase diagram - to show bimodal energy function
%--------------------------------------------------------------------------
x       = -32:1/16:32;
dx      =  x/2 + 25*x./(1 + x.^2);
dxdt    = -x/2 + 16*x./(1 + x.^2) ;
 
subplot(2,1,1)                                 
plot(x,x*0,':',x,dxdt)
axis square
xlabel('state')
ylabel('velocity')
title('phase diagrmas')
 
subplot(2,1,2)
plot(x,dx,x,x,':')
axis square
xlabel('state')
ylabel('image')
 
drawnow
 
 
% temporal correlations
%--------------------------------------------------------------------------
M(1).E.s  = 1/2;
M(1).E.n  = 8;
M(1).E.d  = 3;
M(1).E.nD = 1;
M(1).E.nN = 1;
 
% model specification - 1st level
%--------------------------------------------------------------------------
f      = '(16*x/(1 + x^2) - x/2 + v*2)/8';
g      = '(x^2)/16';
M(1).x = 1;
M(1).f = inline(f,'x','v','P');
M(1).g = inline(g,'x','v','P');
M(1).V = exp(2);
M(1).W = exp(16);
 
% model specification - 2nd level
%--------------------------------------------------------------------------
M(2).v = 0;
M(2).V = 1/8;
 
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
