function spm_ind_demo
% Demo for models of induced frequency responses and nonlinear coupling
%==========================================================================
%
% This demo shows how the nonlinearity in a neural-mass model's sigmoid
% activation function can induce cross-frequency coupling in the outputs.
% In this demo [gamma] frequencies in the response are induced that are
% not in the input. This is the basis of DCM for induced response where
% nonlinear coupling is modelled as coupling between frequency modes. See
% Chen et al for further details:
% 
% Dynamic causal modelling of induced responses
% 
% C.C. Chen, S.J. Kiebel, and K.J. Friston
% 
% ABSTRACT
% 
% This paper describes a dynamic causal model (DCM) for induced or spectral
% responses as measured with the electroencephalogram (EEG) or the
% magnetoencephalogram (MEG). We model the time-varying power, over a range
% of frequencies, as the response of a distributed system of coupled
% electromagnetic sources to a spectral perturbation. The model parameters
% encode the frequency response to exogenous input and coupling among
% sources and different frequencies. The Bayesian inversion of this model,
% given data enables inferences about the parameters of a particular model
% and allows us to compare different models, or hypotheses. One key aspect
% of the model is that it differentiates between linear and nonlinear
% coupling; which correspond to within and between-frequency coupling
% respectively. To establish the face validity of our approach, we generate
% synthetic data and test the identifiability of various parameters to
% ensure they can be estimated accurately, under different levels of noise.
% We then apply our model to EEG data from a face-perception experiment, to
% ask whether there is evidence for nonlinear coupling between early visual
% cortex and fusiform areas.
% 
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_ind_demo.m 5934 2014-03-28 15:03:00Z karl $


% number of regions in coupled map lattice
%--------------------------------------------------------------------------
n     = 1;
 
% specify network (connections)
%--------------------------------------------------------------------------
if n > 1
    A{1} = diag(ones(n - 1,1),-1);
else
    A{1} = 0;
end
A{2}  = A{1}';
A{3}  = sparse(n,n);
B     = {};
C     = sparse(1,1,1,n,1);
 
% create LFP model
%--------------------------------------------------------------------------
M.dipfit.type  = 'LFP';
M.dipfit.model = 'LFP';
M.dipfit.Ns    = n;
M.dipfit.Nc    = n;

% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(A,B,C);
[pE,pC] = spm_L_priors(M.dipfit,pE,pC);

M.f   = 'spm_fx_lfp';
M.g   = 'spm_gx_erp';
M.x   = sparse(n,13);
M.pE  = pE;
M.pC  = pC;
M.m   = size(C,2);
M.n   = n*13;
M.l   = size(pE.L,1);
 
% Integrate system to see response
%--------------------------------------------------------------------------
N     = 256;                          % number of samples
U.dt  = 8/1000;                       % sampling interval
pst   = (1:N)*U.dt;                   % peristimulus time
t     = 64;                           % sample window for WFT
cpt   = 1:1/8:16;                     % cycles per window         
w     = cpt./(t*U.dt);                % Hz
 
% input 
%==========================================================================
U.u   = sparse(64:128,1,1,N,1)*64 + randn(N,1)*4;        % noisy burst
U.u   = sparse(64:128,1,1,N,1).*sin(2*pi*16*pst(:));     % pure Hz - low
U.u   = sparse(64:128,1,1,N,1).*sin(2*pi*16*pst(:))*128; % pure Hz - high
 
% response
%--------------------------------------------------------------------------
fprintf('Generating data; please wait\n')
LFP   = spm_int_ode(pE,M,U);
 
% display
%==========================================================================
spm_figure('GetWin','Figure 1');
 
% input - time
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(pst,U.u)
axis square tight
xlabel('time (s)')
ylabel('activity')
title('input')
 
% Input - time-frequency
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc(pst,w,abs(spm_wft(U.u,cpt,t)));
axis square xy
xlabel('time (s)')
ylabel('frequency')
title('input')
 
% LFP - time
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(pst,LFP)
axis square tight
xlabel('time (s)')
ylabel('activity')
title('response')
 
% LFP - time-frequency
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(pst,w,abs(spm_wft(LFP,cpt,t)));
axis square xy
xlabel('time (s)')
ylabel('frequency')
title('response')
drawnow
