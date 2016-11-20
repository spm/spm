function spm_seizure_demo
% Demo routine for local field potential models
%==========================================================================
% 
% This routine illustrates how one can model induced responses (e.g.,
% seizure onset in terms of exogenously forced changes in model parameters -
% (e.g., recurrent inhibitory connections in a canonical microcircuit
% model. This calls on extra parameters X and Y. X couples input to
% parameters, while Y couples hidden states to parameters.  Here we use
% exogenous input to change the parameters and the ensuing Jacobian to
% elicit fast gamma activity.
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_seizure_demo.m 6937 2016-11-20 12:30:40Z karl $ 
 

% Model specification
%==========================================================================
rng('default')

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
Nc    = 1;
Ns    = 1;
options.spatial  = 'LFP';
options.model    = 'TFM';
options.analysis = 'TFA';

% sub-population showing effect
%--------------------------------------------------------------------------
pop = 2;

M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
M.Hz           = 1:64;

 
% get priors
%--------------------------------------------------------------------------
pE     = spm_dcm_neural_priors({0 0 0},{},0,options.model);
pE     = spm_L_priors(M.dipfit,pE);
pE     = spm_ssr_priors(pE);
x      = spm_dcm_x_neural(pE,options.model);

% eliminate channel noise and make innovations white
%--------------------------------------------------------------------------
pE.a   = [ 0; 0];                  % log amplitude and f^(-a) exponent
pE.b   = [-8; 0];                  % log amplitude and f^(-a) exponent
pE.c   = [-8; 0];                  % log amplitude and f^(-a) exponent


% exogenous input-dependent parameters
%==========================================================================        
np     = length(spm_vec(pE));
nx     = length(spm_vec(x ));
nu     = size(pE.C,2);
i      = spm_fieldindices(pE,'G');
pE.X   = sparse(i(pop),1,1,np,nu);
pE.Y   = sparse(np,nx);
u      = sparse(1,nu);

% create LFP model
%--------------------------------------------------------------------------
M.f    = 'spm_fx_tfm';
M.g    = 'spm_gx_erp';
M.h    = 'spm_fx_cmc_tfm';
M.x    = x;
M.n    = nx;
M.pE   = pE;
M.m    = nu;
M.l    = Nc;
 
% Volterra Kernels and transfer functions
%==========================================================================
spm_figure('GetWin','Volterra kernels and transfer functions');

 
% augment and bi-linearise (with delays)
%--------------------------------------------------------------------------
[f,J,D]       = spm_fx_tfm(x,u,pE,M);
M.u           = sparse(Ns,1);
M.D           = D;
[M0,M1,L1,L2] = spm_bireduce(M,pE);


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


% compute transfer functions for different (self) inhibitory connections
%--------------------------------------------------------------------------
p     = linspace(-2,2,64);
for i = 1:length(p)
    P        = pE;
    P.G(pop) = p(i);
    [G w]    = spm_csd_mtf(P,M);
    GW(:,i)  = abs(G{1});
end

subplot(2,2,3)
plot(w,GW)
xlabel('frequency {Hz}')
title('transfer function','FontSize',16)
axis square

subplot(2,2,4)
imagesc(p,w,log(GW))
title('transfer functions','FontSize',16)
ylabel('Frequency')
xlabel('Inhibitory connection','FontSize',16)
axis xy
axis square

% Integrate system to see response (time-frequency)
%==========================================================================
spm_figure('GetWin','spontaneous fluctuations');


% remove M.u to invoke exogenous inputs
%--------------------------------------------------------------------------
try, M = rmfield(M,'u'); end

N     = 512;
U.dt  = 4/1000;
t     = (1:N)'*U.dt;
U.u   = sparse(N,M.m);

% exogenous input
%--------------------------------------------------------------------------
U.u(:,1) = tanh((t - 1)*8)/2;
M.W      = inv(diag(sparse(1,1,1,1,M.n) + exp(-32)));
LFP      = spm_int_sde(pE,M,U);
 
% plot
%--------------------------------------------------------------------------
subplot(4,1,1)
plot(t,U.u)
xlabel('time (s)')
title('Exogenous input','FontSize',16)
spm_axis tight

% LFP – random fluctuations
%--------------------------------------------------------------------------
subplot(4,1,2)
plot(t,LFP)
xlabel('time (s)')
title('LFP response','FontSize',16)
spm_axis tight
 
% time-frequency
%--------------------------------------------------------------------------
W     = 128;
TFR   = spm_wft(LFP,w*W*U.dt,W);
subplot(4,1,3)
imagesc(t,w,spm_en(abs(TFR)));
title('time-frequency response','FontSize',16)
axis  xy
xlabel('time (s)')
ylabel('Hz')
drawnow

% now integrate a generative model to simulate a time frequency response
%==========================================================================
M.f       = M.h;
M         = rmfield(M,'h');
csd       = spm_csd_int(pE,M,U);

% predicted time frequency response
%--------------------------------------------------------------------------
subplot(4,1,4)
imagesc(t,w,spm_en(abs(csd{1}')));
title('Predicted response','FontSize',16)
axis xy
xlabel('time (s)')
ylabel('Hz')

return

 
 