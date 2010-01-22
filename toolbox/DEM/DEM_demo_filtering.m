% State-space demo routine comparing Bayesian filtering and DEM:  The
% system here is chosen to highlight changes in conditional moments
% (including precision) induced by nonlinearities in the model.  A
% comparative evaluation is provided using extended Kalman filtering and
% particle filtering. Crucially, DEM and particle filtering deal gracefully
% with nonlinearities, in relation to Kalman filtering.

%==========================================================================
clear M
Fgraph    = spm_figure('GetWin','Graphics');
 
% temporal correlations
%--------------------------------------------------------------------------
M(1).E.s  = 1/32;
M(1).E.nD = 4;
M(1).E.linear = 3;
 
% model specification - 1st level
%--------------------------------------------------------------------------
f       = 'exp(v) - P*x';
g       = '(x.^2)/5';
M(1).x  = 1;
M(1).f  = inline(f,'x','v','P');
M(1).g  = inline(g,'x','v','P');
M(1).pE = log(2);
M(1).V  = exp(4);
 
% model specification - 2nd level
%--------------------------------------------------------------------------
M(2).v  = 0;
M(2).V  = 2;
 
% generate data (ouput)
%--------------------------------------------------------------------------
T       = 32;
U       = 1/2 + sin(pi*[1:T]/16);
DEM     = spm_DEM_generate(M,U);
spm_DEM_qU(DEM.pU)
 
% DEM
%--------------------------------------------------------------------------
DEM.U   = U;
DEM     = spm_DEM(DEM);
d_x     = DEM.qU.x{1};
t_x     = DEM.pU.x{1};
 
% overlay true values
%--------------------------------------------------------------------------
spm_DEM_qU(DEM.qU,DEM.pU)
drawnow

% EKF
%==========================================================================
[e_x,Pk] = spm_ekf(DEM.M,DEM.Y);
 
% PF
%==========================================================================
[p_x,Pp] = spm_pf(DEM.M,DEM.Y,U);
 
 
% graphics comparing PF, EKF and DEM
%--------------------------------------------------------------------------
figure(Fgraph)

subplot(2,1,1)
plot([1:T],t_x,'r-',[1:T],d_x,'k-',[1:T],p_x,'k-.',[1:T],e_x,'k:')
legend({'true','DEM','PF','EKF'})
axis square
xlabel('time')
title('Conditional expectation (hidden states)')
 
% graphics conditional covariance
%--------------------------------------------------------------------------
subplot(2,1,2)
Cd     = spm_cat(DEM.qU.S);
Ck     = spm_cat(Pk);
Cp     = spm_cat(Pp);
Ck     = Ck(2:T);
plot(1:T,Cd,'k-',1:T,Cp,'k-.',2:T,Ck,'k:')
legend({'DEM','PF','EKF'})
axis square
xlabel('time')
title('Conditional covariance (hidden states)')
 
 
return
 
 
% repeat for several realizations
%==========================================================================
for i = 1:8
    
    z    = randn(1,T)/sqrt(2);
    DEM  = spm_DEM_generate(M,U + z);
    
    DEM  = spm_DEM(DEM);
    t_x  = DEM.pU.x{1};
    d_x  = DEM.qU.x{1};
    e_x  = spm_ekf(DEM.M,DEM.Y);
    p_x  = spm_pf(DEM.M,DEM.Y,U);
 
    
    SSE(i,1) = sum(sum((e_x - t_x).^2));
    SSE(i,2) = sum(sum((p_x - t_x).^2));
    SSE(i,3) = sum(sum((d_x - t_x).^2));
    
end
 
clf
subplot(2,1,1)
semilogy(SSE','k:'), hold on
semilogy(SSE','k.','Markersize',16), hold off
set(gca,'Xtick',[1 2 3],'XLim',[0 4])
set(gca,'Xticklabel',{'EKF','PF','DEM'})
title('sum of squared error (hidden states)')
axis square
