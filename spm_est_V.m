function [h] = spm_est_V(SPM,c)
% Test routine to evaluate non-sphericity correction (ReML Whitening)
% FORMAT [h] = spm_est_V(SPM,c)
% SPM    - structure containing generic analysis details
% c      - number of contrasts to simulate (default = 4)
%
% h      - hyperparameter estimates
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_est_V.m 5219 2013-01-29 17:07:07Z spm $
 
% get data and model
%==========================================================================
spm_figure('GetWin','Figure');
 
% default number of contrasts
%--------------------------------------------------------------------------
try, c; catch, c = 4; end
 
% check filenames
%--------------------------------------------------------------------------
SPM.xY.VY = spm_check_filename(SPM.xY.VY);
 
 
% get data from significant voxels
%--------------------------------------------------------------------------
N     = 4000;                               % number of voxels
Vspm  = SPM.xCon(1).Vspm;                   % get first SPM
XYZ   = SPM.xVol.XYZ;
F     = spm_sample_vol(Vspm,XYZ(1,:),XYZ(2,:),XYZ(3,:),0);
[F,i] = sort(F,2,'descend');
XYZ   = XYZ(:,i(1:16000));                  % voxels for t-test
rpv   = SPM.xVol.R(end)/SPM.xVol.S;         % resels per voxel
 
 
% get data and covariance
%--------------------------------------------------------------------------
Y     = spm_get_data(SPM.xY.VY,XYZ);
m     = size(Y,1);                          % number of scans
% Y   = spm_null_data(Y,SPM);               % uncomment for null data
 
% Data covariance
%--------------------------------------------------------------------------
C     = cov(Y(:,1:N)');
 
% covariance components (a mixture of exponentials)
%==========================================================================
dt    = SPM.xY.RT;                          % TR (seconds)
T     = (0:(m - 1))*dt;                     % time
d     = 2.^(floor(log2(dt/4)):log2(64));    % time constants (seconds)
QQ    = {};                                 % dictionary of components
for i = 1:length(d)
    for j = 0:1
        QQ{end + 1} = toeplitz((T.^j).*exp(-T/d(i)));
    end
end
 
Q{1}  = QQ(1);                              % white (almost)
Q{2}  = QQ([1,5]);                          % standard (roughly)
Q{3}  = QQ(1:end);                          % full (exactly)
 
 
% estimate serial correlations (and perform null t-tests)
%==========================================================================
 
% get design and augment with drift terms
%--------------------------------------------------------------------------
t     = -6:1/8:6;                           % range of t-values to plot
X     = SPM.xX.X;
try
    X = [X SPM.xX.K.X0];
end
 
% add simulated effects
%--------------------------------------------------------------------------
X     = [spm_conv(randn(size(X,1),c),8/dt,0) X];
 
% Residual forming matrix and scale data covariance
%--------------------------------------------------------------------------
[m,n] = size(X);
R     = speye(m,m) - X*spm_pinv(X);
C     = C*trace(R*R)/trace(R*C*R);
for q = 1:length(Q)
    
    % ReML and whitening matrix (W)
    %----------------------------------------------------------------------
    [V,h]   = spm_reml(C,X,Q{q},1,1,1,0,4);
    W       = spm_inv(spm_sqrtm(V));
    
    % scales for plotting later and effective degrees of freedom
    %----------------------------------------------------------------------
    W       = W*sqrt(trace(R*R)/trace(R*W*C*W*R));
    edf(q)  = trace(R*V)^2/trace(R*V*R*V);
    
    % empirical t-distribution
    %----------------------------------------------------------------------
    for i = 1:c
        [T,tdf] = spm_ancova(W*X,speye(m,m),W*Y,sparse(i,1,1,n,1));
        try
            Tpdf{q} = Tpdf{q} + hist(T,t);
        catch
            Tpdf{q} = hist(T,t);
        end
    end
    
end
 
 
% Fourier transforms
%--------------------------------------------------------------------------
S  = spm_sqrtm(C);
g  = [  sum(abs(fft(full(R)).^2),2)];       % residual forming matrix
g  = [g sum(abs(fft(full(R*S)).^2),2)];     % residuals unwhitened
g  = [g sum(abs(fft(full(R*W*S)).^2),2)] ;  % residuals whitened
 
subplot(2,2,1)
i  = fix(2:m/2);
w  = (1:length(i))/(2*length(i)*dt);
plot(w,g(i,:))
title('Spectral density','FontSize',16)
xlabel('Frequency (Hz)')
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
lag   = -32:32;
i     = lag + fix(m/2);
plot(lag*dt,r(i,:))
title('Auto-covariance function','FontSize',16)
xlabel('lag (seconds)')
ylabel('covariance')
axis square

% return unless contrasts have been simulated
%--------------------------------------------------------------------------
if ~c, return, end
 

% plot FPR above a t-threshold u = 3
%==========================================================================
TPDF  = spm_Tpdf(t,tdf(2));
TPDF  = sum(Tpdf{1})*TPDF/sum(TPDF);
u     = find(abs(t) > 3);
FPR   = sum(TPDF(u));
for q = 1:length(Tpdf)
    fpr(q) = sum(Tpdf{q}(u));
end
 
 
% number of covariance components and effective d.f.
%--------------------------------------------------------------------------
for i = 1:length(Q)
    str{i} = sprintf('%i (%.0f)',length(Q{i}),edf(i));
end
 
% plot in terms of resels (under Poisson assumptions)
%--------------------------------------------------------------------------
subplot(2,2,4)
spm_plot_ci(fpr(:)*rpv,fpr(:)*rpv), hold on
plot([0 length(Tpdf)],[FPR FPR]*rpv,'LineWidth',4), hold off
title('False positive rates u = 3','FontSize',16)
ylabel('Resolution elements')
xlabel('Components (d.f.)')
set(gca,'XTickLabel',str)
axis square
 
% null distributions
%--------------------------------------------------------------------------
subplot(2,2,3)
semilogy(t,Tpdf{end},'r',t,TPDF,'b-.',t,Tpdf{1},'g')
title('t-distributions','FontSize',16)
xlabel('t-value')
ylabel('log-frequency')
axis square
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


