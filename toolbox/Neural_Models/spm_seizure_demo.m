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
% $Id: spm_seizure_demo.m 4718 2012-04-19 15:34:45Z karl $ 
 

% Model specification
%==========================================================================

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'CMC';
options.analysis = 'TFA';
dipfit.model = options.model;
dipfit.type  = options.spatial;
dipfit.Nc    = Nc;
dipfit.Ns    = Ns;

 
% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_dcm_neural_priors({0 0 0},{},[1 0],options.model);
[pE,pC] = spm_L_priors(dipfit,pE,pC);
[pE,pC] = spm_ssr_priors(pE,pC);
[x,f]   = spm_dcm_x_neural(pE,options.model);

% eliminate channel noise and make innovations white
%--------------------------------------------------------------------------
pE.a    = [  0; -16];                  % log amplitude ang f^(-a) exponent
pE.b    = [-32; -32];                  % log amplitude ang f^(-a) exponent
pE.c    = [-32; -32];                  % log amplitude ang f^(-a) exponent


% exogenous input-dependent parameters
%==========================================================================
np      = length(spm_vec(pE));
nx      = length(spm_vec(x ));
nu      = size(pE.C,2);
i       = spm_fieldindices(pE,'G');
j       = 4;
pE.X    = sparse(i(j),2,1,np,nu);
pC.X    = sparse(np,nu);
pE.Y    = sparse(np,nx);
pC.Y    = sparse(np,nx);
u       = sparse(1,nu);

% create LFP model
%--------------------------------------------------------------------------
M.f     = 'spm_fx_tfm';
M.g     = 'spm_gx_erp';
M.h     = f;
M.x     = x;
M.n     = nx;
M.pE    = pE;
M.m     = nu;
M.l     = Nc;
 
% Volterra Kernels and transfer functions
%==========================================================================
spm_figure('GetWin','Volterra kernels and transfer functions');

 
% augment and bi-linearise (with delays)
%--------------------------------------------------------------------------
[f,J,D]       = spm_fx_tfm(x,u,pE,M);
M.u           = sparse(Ns,1);
[M0,M1,L1,L2] = spm_bireduce(M,pE,D);


% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
N          = 64;
dt         = 1/1000;
t          = (1:N)*dt*1000;
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
B     = 2.3;
p     = linspace(-2,B,64);
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
imagesc(p,w,log(GW))
title('transfer functions','FontSize',16)
ylabel('Frequency')
xlabel('Inhibitory connection','FontSize',16)
axis xy


% Integrate system to see response (time-frequency)
%==========================================================================
spm_figure('GetWin','spontaneous fluctuations');


% remove M.u to invoke exogenous inputs
%--------------------------------------------------------------------------
M     = rmfield(M,'u');
N     = 512;
U.dt  = 4/1000;
t     = (1:N)'*U.dt;
U.u   = sparse(N,M.m);

% exogenous input
%--------------------------------------------------------------------------
U.u(:,1) = exp(-(t - 1).^2*16)*2;
U.u(:,2) = tanh((t - 1/2)*4)*B;
X        = spm_int_J(pE,M,U);
M.W      = inv(diag(sparse(1,1,exp(2),1,M.n) + exp(-32)));
LFP      = spm_int_sde(pE,M,U);
 
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
W     = 128;
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

 
 