function spm_lfp_demo
% Demo routine for local field potential models
%==========================================================================
%
% This is a generic demonstration of neural-mass models that illustrates
% various impulse response behaviours. It is meant to show how to specify
% a neural-mass model, examine its response properties using Volterra
% kernels and transfer functions and generate electrophysiological and
% hemodynamic responses from the same model. It is anticipated that people
% will go through the code to see how the routines relate to each other.
%
% This demo contains a linear stability analysis, which can be useful for
% identifying useful domains of parameter space (here the inhibitory time-
% constant)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lfp_demo.m 6937 2016-11-20 12:30:40Z karl $


% Model specification
%==========================================================================

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
n     = 1;

% specify network (connections)
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

M.dipfit.model = 'LFP';
M.dipfit.type  = 'LFP';
M.dipfit.Nc    = n;
M.dipfit.Ns    = n;

% get priors
%--------------------------------------------------------------------------
pE    = spm_dcm_neural_priors(A,B,C,M.dipfit.model);  % neuronal priors
pE    = spm_L_priors(M.dipfit,pE);                    % spatial  priors
pE    = spm_ssr_priors(pE);                           % spectral priors

[x,f] = spm_dcm_x_neural(pE,M.dipfit.model);


% create LFP model
%--------------------------------------------------------------------------
M.f   = f;
M.g   = 'spm_gx_erp';
M.x   = x;
M.pE  = pE;
M.m   = size(C,2);
M.n   = n*13;
M.l   = size(pE.L,1);



% create BOLD model
%--------------------------------------------------------------------------
[hE,hC] = spm_hdm_priors(1,5);
hE(end) = 1;

% model
%--------------------------------------------------------------------------
clear H
H.f     = 'spm_fx_hdm';
H.g     = 'spm_gx_hdm';
H.x     = [0 0 0 0]';
H.pE    = hE;
H.pC    = hC;
H.m     = 1;
H.n     = 4;
H.l     = 1;

% Volterra Kernels
%==========================================================================
spm_figure('GetWin','Volterra kernels');

% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L1,L2] = spm_bireduce(M,M.pE);

% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
N          = 128;
dt         = 1/1000;
t          = (1:N)*dt*1000;
w          = ((1:N) - 1)/N/dt;
[K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);
G1         = fft(K1);
G2         = fft2(K2);
i          = 1:(N/4);

subplot(2,2,1)
plot(t,K1)
title('1st-order Volterra kernel','FontSize',16)
axis square, xlabel('time (ms)')

subplot(2,2,2)
plot(w(i),abs(G1(i)).^2)
title('1st-order transfer function','FontSize',16)
axis square, xlabel('frequency (Hz)')

subplot(2,2,3)
imagesc(t,t,K2)
title('2nd-order Volterra kernel','FontSize',16)
axis square
xlabel('time (ms)')

subplot(2,2,4)
imagesc(w(i),w(i),abs(G2(i,i)).^2)
title('2nd-order transfer function','FontSize',16)
axis square, xlabel('frequency (Hz)')
drawnow


% Integrate system to see response (time-frequency and hemodynamic)
%==========================================================================
spm_figure('GetWin','induced and haemodynamic responses');

N     = 2048;
U.dt  = 8/1000;
U.u   = 32*(sparse(128:512,1,1,N,M.m) + randn(N,M.m)/16);
t     = (1:N)*U.dt;
LFP   = spm_int_L(pE,M,U);

% input
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(t,U.u)
axis square
xlabel('time (s)')
title('Exogenous input')

% LFP
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(t,LFP)
axis square
xlabel('time (s)')
title('LFP response')


% time-frequency
%--------------------------------------------------------------------------
W     = 512;
w     = 4:1/4:32;
cpW   = w*W*U.dt;
subplot(2,2,3)
imagesc(t,w,abs(spm_wft(LFP(:,1),cpW,W)));
title('time-frequency response')
axis square xy
xlabel('time (s)')
ylabel('Hz')

% Use response to drive a hemodynamic model
%--------------------------------------------------------------------------
U.u     = LFP;
BOLD    = spm_int_L(hE,H,U);

subplot(2,2,4)
plot(t,BOLD)
title('BOLD response')
axis square
xlabel('time (s)')
drawnow



% Stability analysis (over excitatory and inhibitory time constants)
%==========================================================================
fprintf('Stability analysis - please wait\n')
P1    = linspace(-4,2,32);
P2    = linspace(-2,1,32);
for i = 1:length(P1)
    for j = 1:length(P2)
        P        = M.pE;
        P.T(:,1) = P.T(:,1) + P1(i);
        P.T(:,2) = P.T(:,2) + P2(j);
        S        = spm_ssm2s(P,M);
        S        = S(abs(imag(S)) > 2*2*pi & abs(imag(S)) < 64*2*pi);
        try
            [s k]   = max(real(S));
            HZ(i,j) = abs(imag(S(k)))/(2*pi);
            LE(i,j) = s;
        end
    end
end

p1  = 4*exp(P1);
p2  = 16*exp(P2);

% graphics
%--------------------------------------------------------------------------
spm_figure('GetWin','bifurcation analysis');

subplot(2,2,1)
surf(p1,p2,LE')
shading interp
axis square
ylabel('inhibitory time constant (ms)')
xlabel('excitatory time constant (ms)')
title('stability','FontSize',16)

subplot(2,2,2)
contour(p2,p1,LE)
axis square
xlabel('inhibitory time constant')
ylabel('excitatory time constant')
title('stability','FontSize',16)

subplot(2,2,3)
surf(p1,p2,HZ')
shading interp
axis square
ylabel('inhibitory time constant')
xlabel('excitatory time constant')
title('Frequency','FontSize',16)

subplot(2,2,4)
contour(p2,p1,full(HZ),[8 12 16 30],'k');
axis square
xlabel('inhibitory time constant')
ylabel('excitatory time constant')
title('Frequency','FontSize',16)
drawnow



% transfer functions
%==========================================================================
spm_figure('GetWin','transfer functions');


% compute transfer function
%--------------------------------------------------------------------------
pE    = M.pE;
[G w] = spm_lfp_mtf(pE,M);

subplot(2,1,1)
plot(w,G{1})
axis square
xlabel('frequency {Hz}')
title('transfer function')
drawnow


% compute transfer functions for different inhibitory time constants
%--------------------------------------------------------------------------
p     = log((1:64)/32);
for i = 1:length(p)
    pE.T(2) = p(i);
    G       = spm_lfp_mtf(pE,M);
    GW(:,i) = G{1};
end

subplot(2,2,3)
imagesc(16*exp(p),w,GW)
ylabel('Frequency')
xlabel('Inhibitory time constant (ms)')

subplot(2,2,4)
plot(w,GW)
xlabel('Frequency')
ylabel('g(w)')
drawnow


% Integrate system to see Transient response (with noise)
%==========================================================================
spm_figure('GetWin','Time frequency analysis');

pE.T(2) = log(2);
N     = 1024;
U.dt  = 1/1000;
U.u   = 8*(exp(-((1:N)' - N/4).^2/(2*32^2)) + randn(N,1)/4);
t     = [1:N]*U.dt;
LFP   = spm_int_L(pE,M,U);

% LFP
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(t*1000,U.u)
title('input')
axis square
xlabel('time (ms)')

% LFP
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(t*1000,LFP)
title('response')
axis square
xlabel('time (ms)')

% time-frequency
%--------------------------------------------------------------------------
W     = 512;
w     = 4:1/4:32;
cpW   = w*W*U.dt;
subplot(2,2,3)
imagesc(t*1000,w,abs(spm_wft(LFP(:,1),cpW,W)));
title('time-frequency response')
axis square xy
xlabel('time (ms)')
