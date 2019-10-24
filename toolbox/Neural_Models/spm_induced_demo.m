% function spm_induced_demo
% Demo routine for induced responses
%==========================================================================
%
% This demonstration illustrates the generative or forward model used for
% time frequency responses - in other words, a biophysically plausible
% dynamic causal model for induced responses. The basic idea is to
% integrate a neural mass model to obtain expected hidden neuronal
% states produced by an unknown (parameterised) exogenous input. The states 
% are used as the expansion point for a linear perturbation analysis of
% the frequency response properties that are local in peristimulus time.
% The ensuing spectra (induced complex cross spectra) are then convolved
% with a wavelet window to generate predictions of a conventional time
% frequency (wavelet) transform. Crucially, these predictions are complex
% and can be used to characterise delays - in terms of cross covariance
% functions. Nonlinearities in the neural mass model mean that the spectral
% responses caused by random neuronal fluctuations are state dependent and
% therefore change with the expected hidden states over peristimulus time.
%
% This routine first creates a simple - two source - generative model using
% a canonical microcircuit architecture and convolution based dynamics. It
% then produces predictions of induced responses to a short and sustained
% input to the first source - as measured by two local field potential
% recordings at each source. Exactly the same model is then integrated in
% time,  using (serially correlated) random fluctuations to drive each
% source (in addition to the exogenous input). This is repeated over 16
% trials and the simulated responses are characterised in terms of a
% wavelet transform - to produce complex cross spectral  data features. 
% These are shown graphically with their analytic predictions from the 
% generative model.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_induced_demo.m 7679 2019-10-24 15:54:07Z spm $
 
 
% Model specification
%==========================================================================
rng('default')
clear spm_fp_cmc_tfm
 
% number of regions in coupled map lattice
%--------------------------------------------------------------------------
Nc    = 2;                                       % number of channels
Ns    = 2;                                       % number of sources
options.spatial  = 'LFP';
options.model    = 'TFM';
options.analysis = 'TFM';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
 
% extrinsic connections (forward an backward)
%--------------------------------------------------------------------------
A{1} = [0 0; 1 0];
A{2} = [0 1; 0 0];
A{3} = [0 0; 0 0];
C    = [1 0;0 1];

% get priors
%--------------------------------------------------------------------------
pE    = spm_dcm_neural_priors(A,{},C,options.model);
pE    = spm_L_priors(M.dipfit,pE);
pE    = spm_ssr_priors(pE);
x     = spm_dcm_x_neural(pE,options.model);

% extrinsic connections
%--------------------------------------------------------------------------
pE.A{1}   = pE.A{1} + 1;
pE.A{2}   = pE.A{2} + 1;
pE.A{3}   = pE.A{3} - 1;
pE.A{4}   = pE.A{4} - 1;

% supress noise
%--------------------------------------------------------------------------
pE.a(1,:) = -2;
pE.b(1,:) = -4;
pE.c(1,:) = -4;
 
% orders and model
%==========================================================================
nx    = length(spm_vec(x));
nu    = size(pE.C,2);


% create forward model
%--------------------------------------------------------------------------
M.g   = 'spm_gx_erp';
M.f   = 'spm_fx_cmc_tfm';
M.h   = 'spm_fp_cmc_tfm';
M.x   = x;
M.n   = nx;
M.pE  = pE;
M.m   = nu;
M.l   = Nc;
M.Hz  = 2:64;
M.Rft = 4;
M.ds  = 4;
 
% solve for steady state
%--------------------------------------------------------------------------
M.x   = spm_dcm_neural_x(pE,M);
 
% Integrate system to see response (time-frequency)
%==========================================================================
 
% invoke exogenous inputs
%--------------------------------------------------------------------------
N     = 128;
U.dt  = 6/1024;
pst   = ((1:N)' - N/4)*U.dt;
U.u   = sparse(N,M.m);
 

% exogenous input - a sustained input of about 128 seconds
%--------------------------------------------------------------------------
M.ons    = 64;
U.u(:,1) = spm_erp_u(pst,pE,M);

 
% integrate generative model to simulate a time frequency response
%--------------------------------------------------------------------------
[csd,erp,CSD,mtf,w,pst,x,dP] = spm_csd_int(pE,M,U);


% plot expected responses
%==========================================================================
spm_figure('GetWin','Simulated time-frequency responses');
pst = 1000*pst;
 
subplot(4,2,1)
plot(pst,U.u)
xlabel('peristimulus time (ms)')
title('Exogenous input','FontSize',16)
spm_axis tight
 
% LFP - expectation
%--------------------------------------------------------------------------
subplot(4,2,2)
j   = find(kron(sparse(1,[3 1 5],1,1,8),ones(Ns,1)));
plot(pst,x(j,:))
xlabel('peristimulus time (ms)')
title('Hidden neuronal states','FontSize',16)
spm_axis tight

subplot(4,2,3)
plot(pst,erp{1})
xlabel('peristimulus time (ms)')
title('Evoked response','FontSize',16)
spm_axis tight
 
% LFP - expectation
%--------------------------------------------------------------------------
subplot(4,2,4)
plot(pst,dP{1})
xlabel('peristimulus time (ms)')
title('intrinsic connectivity','FontSize',16)
spm_axis tight


% expected time frequency response (coherence and cross-covariance)
%--------------------------------------------------------------------------
spm_dcm_tfm_image(CSD{1},pst,w,1)


% expected time frequency response
%--------------------------------------------------------------------------
spm_figure('GetWin','transfer functions');

spm_dcm_tfm_transfer(mtf{1},pst,w)


% simulated responses
%==========================================================================
spm_figure('GetWin','Predicted responses');
 
% time-frequency
%--------------------------------------------------------------------------
xY.erp = erp;
xY.csd = csd;
spm_dcm_tfm_response(xY,pst,w)



% Integrate system to simulate responses
%==========================================================================
spm_figure('GetWin','Simulated trials');


% get spectral density of neuronal fluctuations
%--------------------------------------------------------------------------
Hz    = 2:128;
Gu    = spm_csd_mtf_gu(pE,Hz);
 
% simulate Nt trials
%--------------------------------------------------------------------------
M.ds  = 0;                                       % switch to ERP generation
Nt    = 8;
V     = U;
for j = 1:Nt
    
    fprintf('\nsimulating trial %i',j)
    y    = 32;
    while max(abs(spm_vec(y))) > 1
        V.u  = U.u + spm_rand_power_law(Gu,Hz,U.dt,N)*2;
        y    = spm_csd_int(pE,M,V);
    end
    D(:,:,j) = full(real(y{1}));
    
    % LFP - random fluctuations
    %----------------------------------------------------------------------
    subplot(2,2,1)
    plot(pst,D(:,:,j)), hold on
    xlabel('time (s)')
    title('simulated response and ERP','FontSize',16)
    spm_axis tight, axis square; drawnow
    
end
hold off

% LFP - expectations
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(pst*1000,erp{1},'-.'), hold on
plot(pst*1000,mean(D,3)), hold off
xlabel('time (s)')
title('LFP response - expectation','FontSize',16)
spm_axis tight, axis square


% Time frequency response
%--------------------------------------------------------------------------
Nf    = length(M.Hz);
P     = zeros(N,Nf,Nc,Nc);
Q     = zeros(N,Nf,Nc,Nc);
E     = mean(D,3);
for k = 1:Nt
    
    d     = full(double(D(:,:,k)));
    G     = spm_morlet(d - E,w*U.dt,M.Rft);
    for i = 1:Nc
        for j = 1:Nc
            P(:,:,i,j) = (G(:,:,i).*conj(G(:,:,j)));
        end
    end
    Q = Q + P;
end

% scale induced responsest sample variance
%--------------------------------------------------------------------------
Vm    = mean(mean((var(D))));
Vs    = mean(diag(squeeze(mean(squeeze(mean(Q))))));
Q     = Vm*Q/Vs;

% time-frequency
%--------------------------------------------------------------------------
spm_dcm_tfm_image(Q,pst,w,1)

% simulated responses
%==========================================================================
spm_figure('GetWin','Predicted responses - simulated');
 
% time-frequency
%--------------------------------------------------------------------------
xY.erp{1} = E;
xY.csd{1} = Q;
spm_dcm_tfm_response(xY,pst,w)


% compare expected and simulated responses
%==========================================================================
spm_figure('GetWin','Expected and simulated responses');

spm_dcm_tfm_image(csd{1},pst,w,0)
spm_dcm_tfm_image(Q,pst,w,1)
 
 
