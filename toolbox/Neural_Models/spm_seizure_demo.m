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
% $Id: spm_seizure_demo.m 4713 2012-04-10 13:25:39Z karl $ 
 

% Model specification
%==========================================================================

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
type  = 'LFP';
model = 'TFM';
dipfit.model = model;
dipfit.type  = type;
dipfit.Nc    = Nc;
dipfit.Ns    = Ns;

 
% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_dcm_neural_priors({0 0 0},{},[1 1],model);
[pE,pC] = spm_L_priors(dipfit,pE,pC);
[pE,pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,model);

% eliminate channel noise and make innovations might
%--------------------------------------------------------------------------
pE.a    = [  0; -16];
pE.b    = [-16; -16];
pE.c    = [-16; -16];

% create LFP model
%--------------------------------------------------------------------------
M.g     = 'spm_gx_erp';
M.f     = f;
M.x     = x;
M.n     = length(spm_vec(x));
M.pE    = pE;
M.m     = size(pE.C,2);
M.l     = Nc;
 
% Volterra Kernels and transfer functions
%==========================================================================
spm_figure('GetWin','Volterra kernels and transfer functions');
 
% augment and bi-linearise
%--------------------------------------------------------------------------
M.u           = sparse(Ns,1);
[M0,M1,L1,L2] = spm_bireduce(M,M.pE);

% remove M.u to invoke exogenous inputs
%--------------------------------------------------------------------------
M              = rmfield(M,'u');
 
% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
N          = 64;
dt         = 1/1000;
t          = [1:N]*dt*1000;
[K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);
 
subplot(2,2,1)
plot(t,K1(:,:,1))
title('1st-order Volterra kernel','FontSize',16)
axis square
xlabel('time (ms)')
 
subplot(2,2,2)
imagesc(t,t,K2(1:64,1:64,1,1,1))
title('2nd-order Volterra kernel','FontSize',16)
axis square
xlabel('time (ms)')


% compute transfer functions for different inhibitory connections
%--------------------------------------------------------------------------
p     = linspace(-4,1.7,32);
for i = 1:length(p)
    P       = pE;
    P.G(4)  = p(i);
    [G w]   = spm_csd_mtf(P,M);
    GW(:,i) = G{1};
end

subplot(2,2,3)
plot(w,GW)
xlabel('frequency {Hz}')
title('transfer function','FontSize',16)
drawnow


subplot(2,2,4)
imagesc(p,w,GW)
title('transfer functions','FontSize',16)
ylabel('Frequency')
xlabel('Inhibitory connection','FontSize',16)
axis xy


% Integrate system to see response (time-frequency)
%==========================================================================
spm_figure('GetWin','spontaneous fluctuations');

N     = 512;
U.dt  = 8/1000;
t     = (1:N)'*U.dt;
U.u   = sparse(N,M.m);


% input
%--------------------------------------------------------------------------
U.u(:,1) = 1;
U.u(:,2) = tanh((t - 2)*2)*1.7;
X        = spm_int_L(pE,M,U);

R        = U;
R.u(:,1) = R.u(:,1) + spm_conv(randn(N,1)*exp(4),2);
LFP      = spm_int_L(pE,M,R);
 
% plot
%--------------------------------------------------------------------------
subplot(4,1,1)
plot(t,U.u)
xlabel('time (s)')
title('Exogenous input','FontSize',16)
spm_axis tight
 
% LFP – expectation
%--------------------------------------------------------------------------
subplot(4,1,2)
plot(t,X)
xlabel('time (s)')
title('LFP response – expectation','FontSize',16)
spm_axis tight

% LFP – random fluctuations
%--------------------------------------------------------------------------
subplot(4,1,3)
plot(t,LFP)
xlabel('time (s)')
title('LFP response','FontSize',16)
spm_axis tight
 
% time-frequency
%--------------------------------------------------------------------------
W     = 98;
TFR   = spm_wft(LFP,w*W*U.dt,W);
subplot(4,1,4)
imagesc(t,w,abs(TFR));
title('time-frequency response','FontSize',16)
axis  xy
xlabel('time (s)')
ylabel('Hz')
drawnow



% now integrate a generative model to simulate a time frequency response
%==========================================================================
[y,w,t,x] = spm_csd_tfm(pE,M,U);

% plot
%--------------------------------------------------------------------------
spm_figure('GetWin','Simulated time frequency responses');

subplot(4,1,1)
plot(t,U.u)
xlabel('time (s)')
title('Exogenous input','FontSize',16)
spm_axis tight
 
% LFP – expectation
%--------------------------------------------------------------------------
subplot(4,1,2)
plot(t,x)
xlabel('time (s)')
title('Hidden neuronal states','FontSize',16)
spm_axis tight

% predicted time frequency response
%--------------------------------------------------------------------------
subplot(4,1,3)
imagesc(t,w,abs(y{1}'));
title('Time-frequency response','FontSize',16)
axis  xy
xlabel('time (s)')
ylabel('Hz')

return

 
 