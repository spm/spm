function DEM_demo_GF_and_KF
% A demonstration of generalised and Kalman filtering where the number
% of hidden states exceeds the number of variables observed. The metrics of
% performance are the mean sum of squared error and the SSE normalized by
% the posterior precision (NESS). The results of a single time series
% analysis are shown first and then the simulations are repeated under
% linear and nonlinear observation models to compare the relative
% performance of DEM and EKF. The superiority of DEM (generalised filtering)
% over Kalman filtering rests on the optimisation of K - the rate of
% generalised descent on free energy (see code after 'return').
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: DEM_demo_GF_and_KF.m 4804 2012-07-26 13:14:18Z karl $
 
 
% linear generative model
%==========================================================================
 
% The problem: A simple linear convolution model with parameters
%--------------------------------------------------------------------------
n       = 8;                  % number of hidden states
m       = 4;                  % number of observations
P.f     = full(toeplitz(sparse([1 2],1,[-2 1],n,1)));
P.g     = spm_dctmtx(m,n);    % Jacobian and observer matrices
 
              
% preliminaries
%--------------------------------------------------------------------------
N        = 128;               % length of sequence
t        = 1:N;               % time
M(1).E.s = 1;                 % smoothness of fluctuations
M(1).E.n = 4;                 % order of generalised coordinates
M(1).E.K = exp(-1);           % rate of generalized gradient descent
M(1).E.linear = 1; 
 
% Level 1: 
%--------------------------------------------------------------------------
M(1).f   = inline('P.f*x','x','v','P');
M(1).g   = inline('P.g*x','x','v','P');
M(1).x   = zeros(n,1);
M(1).pE  = P;
M(1).V   = exp(4);            % precision of observation noise
M(1).W   = exp(2);            % precision of fluctuations on hidden states
 
   
% create data
%==========================================================================
DEM      = spm_DEM_generate(M,N);
 
% Comparative inversions (variants of generalised filtering)
%==========================================================================
 
% DEM and EKF
%--------------------------------------------------------------------------
DEM     = spm_DEM(DEM);
[EKF S] = spm_ekf(DEM.M,DEM.Y);
  
 
% show prediction errors
%==========================================================================
spm_figure('Getwin','Figure 1');
 
MSE(1,:) = mean((DEM.pU.x{1} - DEM.qU.x{1}).^2);
MSE(2,:) = mean((DEM.pU.x{1} - EKF).^2);
 
for i = 1:N
    E         = DEM.pU.x{1}(:,i) - DEM.qU.x{1}(:,i);
    NEES(1,i) = E'*spm_inv(DEM.qU.S{i})*E;
    E         = DEM.pU.x{1}(:,i) - EKF(:,i);
    NEES(2,i) = E'*spm_inv(S{i})*E;
end
 

% plot MSE
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(t,MSE)
title('MSE','Fontsize',16)
xlabel('time','Fontsize',12)
ylabel('mean square error','Fontsize',12)
axis square
legend('DEM','EKF')
 
% plot NEES
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(t,NEES)
title('NEES','Fontsize',16)
xlabel('time','Fontsize',12)
ylabel('normalised error','Fontsize',12)
axis square
legend('DEM','EKF')
drawnow
 
% plot trajectories
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 2');
subplot(2,1,1)
plot(t,DEM.pU.x{1}(1,:),t,DEM.qU.x{1}(1,:),t,EKF(1,:))
title('First hidden state','Fontsize',16)
xlabel('time','Fontsize',12)
axis square
legend('true','DEM','EKF')
drawnow
 
 
% repeat several times with and without a nonlinear observer
%==========================================================================
N     = 64;                         % length of sequence
for i = 1:4
    
    % DEM and EKF (linear)
    %----------------------------------------------------------------------
    M(1).E.linear = 1;
    M(1).g  = inline('P.g*x','x','v','P');
    DEM     = spm_DEM_generate(M,N);
    DEM     = spm_DEM(DEM);
    EKF     = spm_ekf(DEM.M,DEM.Y);
    
    % mean square error
    %----------------------------------------------------------------------
    mse(i,1) = mean(mean((DEM.pU.x{1} - DEM.qU.x{1}).^2));
    mse(i,2) = mean(mean((DEM.pU.x{1} - EKF).^2));
    
    % DEM and EKF (nonlinear)
    %----------------------------------------------------------------------
    M(1).E.linear = 0;
    M(1).g  = inline('tanh(P.g*x)','x','v','P');
    DEM     = spm_DEM_generate(M,N);
    DEM     = spm_DEM(DEM);
    EKF     = spm_ekf(DEM.M,DEM.Y);
    
    % mean square error
    %----------------------------------------------------------------------
    mse(i,3) = mean(mean((DEM.pU.x{1} - DEM.qU.x{1}).^2));
    mse(i,4) = mean(mean((DEM.pU.x{1} - EKF).^2));
        
end
 

spm_figure('Getwin','DEM');
spm_DEM_qU(DEM.qU,DEM.pU)


% show prediction errors
%==========================================================================
spm_figure('Getwin','Figure 1');
 
% plot trajectories
%--------------------------------------------------------------------------
subplot(2,2,3)
plot(1:N,DEM.pU.x{1}(1,:),1:N,DEM.qU.x{1}(1,:),1:N,EKF(1,:))
title('First hidden state','Fontsize',16)
xlabel('time','Fontsize',12)
axis square
legend('true','DEM','EKF')
 
 
% plot MSE
%--------------------------------------------------------------------------
subplot(2,2,4)
plot([1 2],mse(:,[1 2]),'k.:','MarkerSize',16), hold on
plot([3 4],mse(:,[3 4]),'r.:','MarkerSize',16), hold off
title('linear (black) and nonlinear (red)','Fontsize',16)
set(gca,'XTick',[1 2 3 4],'Xlim',[0 5],'XTickLabel',{'DEM','EKF','DEMn','EKFn'})
ylabel('mean square error','Fontsize',12)
axis square


return

% prediction errors as a function of K
%==========================================================================
M(1).E.linear = 1;
M(1).g        = inline('P.g*x','x','v','P');
DEM           = spm_DEM_generate(M,N);

K     = -8:1/2:4;
for i = 1:length(K)
    
    % DEM and EKF (linear)
    %----------------------------------------------------------------------
    DEM.M(1).E.K = exp(K(i));
    DEM          = spm_DEM(DEM);
    
    % mean square error
    %----------------------------------------------------------------------
    kse(i,1) = mean(mean((DEM.pU.x{1} - DEM.qU.x{1}).^2));
    
end

% show prediction errors
%==========================================================================
spm_DEM_qU(DEM.qU,DEM.pU)
 
% plot trajectories
%--------------------------------------------------------------------------
spm_figure('Getwin','Figure 2');
subplot(2,1,2)
plot(K,log(kse))
title('log(MSE)','Fontsize',16)
xlabel('log(K)','Fontsize',12)
axis square

