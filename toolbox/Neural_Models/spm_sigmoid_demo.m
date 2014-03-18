function spm_sigmoid_demo
% Demo routine for neural mass models and the activation function
% =========================================================================
%
% This demo looks at the role of the sigmoid activation function in shaping
% the impulse response of a neural-mass model. It uses Volterra kernels and
% transfer functions and their dependency on the slope parameter of the
% activation function; It is based on the paper by Marreiros et al :
%
% Population dynamics: variance and the sigmoid activation function
%
% André C. Marreiros, Jean Daunizeau, Stefan J. Kiebel, Karl J. Friston
%
% Wellcome Trust Centre for Neuroimaging, University College London, United
% Kingdom
%
% Abstract
%
% This paper demonstrates how the sigmoid activation function in
% neural-mass models can be understood in terms of the variance or
% dispersion of neuronal states.  We use this relationship to estimate the
% probability density on hidden neuronal states, using non-invasive
% electrophysiological (EEG) measures and dynamic casual modelling.  The
% importance of implicit variance in neuronal states for neural-mass models
% of cortical dynamics is illustrated using both synthetic data and real
% EEG measurement of sensory evoked responses.
%
%________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_sigmoid_demo.m 5922 2014-03-18 20:10:17Z karl $

% relating R amd the variance
%==========================================================================
r     = (1:32)/8;
w     = 0;
X     = -64:1/32:64;
for i = 1:length(r)
    dSdx  = (r(i)*exp(-r(i)*(X - w))./(1 + exp(-r(i)*(X - w))).^2);
    dSdx  = dSdx/sum(dSdx);
    V(i)  = sum(X.^2.*dSdx);
end

spm_figure('GetWin','Figure 1'); clf

subplot(2,1,1)
plot(r,sqrt(V))
xlabel('slope parameter')
ylabel('standard deviation')
title('Slope and dispersion','FontSize',16)
axis square
drawnow


% Kernels and transfer functions
%==========================================================================

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
n     = 1;

% specifc network (connections)
%--------------------------------------------------------------------------
if n > 1
    A{1}  = diag(ones(n - 1,1),-1);
else
    A{1} = 0;
end
A{2}  = A{1}';
A{3}  = sparse(n,n);
B     = {};
C     = sparse(1,1,1,n,1);

% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(A,B,C);
[pE,pC] = spm_L_priors(n,pE,pC);
[pE,pC] = spm_ssr_priors(pE,pC);

% create LFP model
%--------------------------------------------------------------------------
M.dipfit.type = 'LFP';

M.f   = 'spm_fx_lfp';
M.g   = 'spm_gx_erp';
M.x   = sparse(n,13);
M.pE  = pE;
M.pC  = pC;
M.m   = size(C,2);
M.n   = n*13;
M.l   = size(pE.L,1);

% Volterra Kernels
%==========================================================================
spm_figure('GetWin','Figure 2'); clf


% kernels
%--------------------------------------------------------------------------
r = [.8 1.6];
m = length(r);
for i = 1:m
    
    pE            = M.pE;
    pE.R(1)       = log(r(i));
    [M0,M1,L1,L2] = spm_bireduce(M,pE);
    
    % compute kernels (over 128ms)
    %----------------------------------------------------------------------
    N          = 128;
    dt         = 1/1000;
    t          = (1:N)*dt*1000;
    [K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);
    
    subplot(2,m,(i - 1) + 1)
    plot(t,K1)
    title(sprintf('1st-order Volterra kernel: slope = %.1f',r(i)),'FontSize',16)
    axis square
    set(gca,'XLim',[t(1) t(end)])
    xlabel('time (ms)')
    
    subplot(2,m,(i - 1) + m + 1)
    imagesc(t,t,K2(1:N,1:N,1))
    title(sprintf('2nd-order Volterra kernel: slope = %.1f',r(i)),'FontSize',16)
    axis square
    xlabel('time (ms)')
    
end

% transfer functions
%==========================================================================
spm_figure('GetWin','Figure 3'); clf

% compute transfer functions for different slope parameters
%--------------------------------------------------------------------------
clear GW
r     = (1:32)/16;
for i = 1:length(r)
    pE.R(1) = log(r(i));
    [G w]   = spm_lfp_mtf(pE,M);
    GW(:,i) = G{1};
end

subplot(2,1,1)
imagesc(w,r,GW')
xlabel('Frequency')
ylabel('slope parameter')
axis xy square
set(gca,'XLim',[w(1) w(end)])

subplot(2,1,2)
plot(w,GW,':r')
xlabel('Frequency')
ylabel('response')
axis square
drawnow


% Integrate system to see transient
%==========================================================================
spm_figure('GetWin','Figure 4'); clf

N     = 256;
U.dt  = 1/1000;
U.u   = 32*(exp(-((1:N)' - N/8).^2/(2*1^2)) + rand(N,1)/8);
t     = (1:N)*U.dt;
r     = [.8 1.6];
m     = length(r);
for i = 1:m
    
    pE.R(1) = log(r(i));
    LFP     = spm_int_L(pE,M,U);
    
    % LFP
    %----------------------------------------------------------------------
    subplot(2,m,(i - 1) + 1)
    plot(t*1000,LFP)
    title(sprintf('depolarization: slope = %.1f',r(i)),'FontSize',16)
    axis square
    xlabel('time (ms)')
    
    % time-frequency
    %----------------------------------------------------------------------
    W     = 256;
    w     = 4:1/4:64;
    cpW   = w*W*U.dt;
    
    subplot(2,m,(i - 1) + m + 1)
    imagesc(t*1000,w,abs(spm_wft(LFP(:,1),cpW,W)).^2);
    title(sprintf('time-frequency response: slope = %.1f',r(i)),'FontSize',16)
    axis square xy
    xlabel('time (ms)')
    ylabel('frequency')
    
end
