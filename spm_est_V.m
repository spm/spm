function [h] = spm_est_V(SPM)
% Test routine to evaluate non-sphericity correction (ReML Whitening)
% FORMAT spm_est_V(SPM)
 
% SPM    - structure containing generic analysis details
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_est_V.m 5033 2012-11-02 20:59:54Z karl $
 
% get data and model
%==========================================================================
spm_figure('GetWin','Figure');
 
 
% check filenames
%--------------------------------------------------------------------------
SPM.xY.VY = spm_check_filename(SPM.xY.VY);
 
 
% get data from responsive voxels
%--------------------------------------------------------------------------
N     = 4000;                                  % number of voxels
Vspm  = SPM.xCon(1).Vspm;                      % get first SPM
XYZ   = SPM.xVol.XYZ;
F     = spm_sample_vol(Vspm,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
[F,i] = sort(F,2,'descend');
XYZ   = XYZ(:,i(1:16000));                     % voxels for t-test
rpv   = SPM.xVol.R(end)/SPM.xVol.S;            % resels per voxel
 
 
% get data and covariance
%--------------------------------------------------------------------------
Y     = spm_get_data(SPM.xY.VY,XYZ);
m     = size(Y,1);                             % number of scans
% Y   = spm_null_data(Y,SPM);                  % uncomment for null data
 
% Data covariance
%--------------------------------------------------------------------------
C     = cov(Y(:,1:N)');
 
% covariance components (a mixture of exponential)
%==========================================================================
dt    = SPM.xY.RT;
T     = (0:(m - 1))*dt;
a     = [1/32 1/2 -1/2 (1 - 1/8) (1 - 1/32)];
a     = [1/8 1/4 1/2 1 2 4 8 16 32 64];
QQ    = {};
for i = 1:length(a)
    QQ{end + 1} = toeplitz(exp(-T/a(i)));
    QQ{end + 1} = toeplitz(T.*exp(-T/a(i)));
end
 
Q{1} = QQ(1);                                % white (almost)
Q{2} = QQ(1:3);                              % standard (almost)
Q{3} = QQ(1:end);                            % full
 
 
% number of covariance components
%--------------------------------------------------------------------------
for i = 1:length(Q)
    nQ(i) = length(Q{i});
end
 
% estimate serial correlations and perform null t-tests
%==========================================================================
 
 
% get design and augment with drift terms
%--------------------------------------------------------------------------
t     = -6:1/8:6;                            % range of t-values to plot
X     = SPM.xX.X;                   
try
    X = [X SPM.xX.K.X0];
end
 
% add simulated effects
%--------------------------------------------------------------------------
E     = spm_conv(randn(size(X,1),1),8,0);
X     = [E X];
 
% Residual forming matrix and scale data covariance
%--------------------------------------------------------------------------
[m,n] = size(X);
R     = speye(m,m) - X*spm_pinv(X);
C     = C*trace(R*R)/trace(R*C*R);
for q = 1:length(Q)
    
    % ReML and whitening matrix (W)
    %----------------------------------------------------------------------
    [V h]   = spm_reml(C,X,Q{q},1,0,2,0,1/2); 
    W       = spm_inv(spm_sqrtm(V));
    W       = W*sqrt(trace(R*R)/trace(R*W*C*W*R));
    
    DF(q)   = trace(R*V)^2/trace(R*V*R*V);
    
    % empirical t-distribution (removing any global bias)
    %----------------------------------------------------------------------
    [T,df]  = spm_ancova(W*X,speye(m,m),W*Y,sparse(1,1,1,n,1));
    T       = T - mean(T);
    Tpdf{q} = hist(T,t);
    
end
 
 
% Fourier transforms
%--------------------------------------------------------------------------
S = spm_sqrtm(C);
g = [  sum(abs(fft(full(R)).^2),2)];
g = [g sum(abs(fft(full(R*S)).^2),2)];
g = [g sum(abs(fft(full(R*W*S)).^2),2)];
 
subplot(2,2,1)
i    = fix(2:m/2);
w    = (1:length(i))/2;
plot(w,g(i,:))
title('Spectral density','FontSize',16)
xlabel('Frequency cycles per session')
ylabel('power')
axis square
legend({'ideal','unwhitened','whitened'})
 
% correlation functions
%--------------------------------------------------------------------------
for i = 1:size(g,2);
    f      = ifft(g(:,i));
    r(:,i) = real(fftshift(f));
end
subplot(2,2,2)
i  = (-32:32) + fix(m/2);
plot(r(i,:))
title('Auto-covariance function','FontSize',16)
xlabel('lag (TR)')
ylabel('covariance')
axis square
 
 
disp(DF);
 
 
% plot FPR above a t-threshold u = 3
%--------------------------------------------------------------------------
subplot(2,2,4)
TPDF  = spm_Tpdf(t,df(2));
TPDF  = sum(Tpdf{1})*TPDF/sum(TPDF);
u     = find(abs(t) > 3);
FPR   = sum(TPDF(u));
for q = 1:length(Tpdf)
    fpr(q)  = sum(Tpdf{q}(u));
end
 
 
% plot in terms of resels (under Poisson assumptions)
%--------------------------------------------------------------------------
spm_plot_ci(fpr(:)*rpv,fpr(:)*rpv), hold on
plot([0 length(Tpdf)],[FPR FPR]*rpv,'LineWidth',4), hold off
title('False positive rates u = 3','FontSize',16)
ylabel('Resolution elements')
xlabel('Number of covariance components')
set(gca,'XTickLabel',nQ)
axis square
 
 
% null distributions
%--------------------------------------------------------------------------
subplot(2,2,3)
semilogy(t,Tpdf{end},'k',t,TPDF,'k-.',t,Tpdf{1},'r')
title('t-distributions','FontSize',16)
xlabel('t-value')
ylabel('log-frequnecy')
axis square
legend({'whitened - empirical','true','unwhitened - empirical'})
drawnow
 
 
% Null data
%==========================================================================
function Y = spm_null_data(Y,SPM)
 
% get design and augment with drift terms
%--------------------------------------------------------------------------
X     = SPM.xX.X; try, X = [X SPM.xX.K.X0]; end
 
 
% reconstitute with phase-shuffled noise
%--------------------------------------------------------------------------
sig   = X*spm_pinv(X)*Y;
res   = Y - sig;
res   = spm_phase_shuffle(res);
res   = spm_conv(randn(size(res)),2)*std(res(:));
Y     = sig + res;
