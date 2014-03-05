% function spm_dcm_Granger_Jacobian
% Demo routine for induced responses
%==========================================================================
%
% This routine illustrates the derivation of spectral Granger causal
% measures from the inversion of a simple state-space DCM paramterised in
% explcicity in terms of its Jacobian (null model).
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_Granger_Jacobian.m 5908 2014-03-05 20:31:57Z karl $
 
 
% Process (model) specification
%==========================================================================
rng('default')
 
% number of regions
%--------------------------------------------------------------------------
Nc    = 2;                                       % number of channels
Ns    = 2;                                       % number of sources
ns    = 2*128;                                   % sampling frequency
dt    = 1/ns;                                    % time bins
Hz    = 1:128;                                   % frequency
p     = 16;                                      % autoregression order
options.spatial  = 'LFP';
options.model    = 'CMC';
options.analysis = 'CSD';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;
M.pF.D         = [1 4];
 
% extrinsic connections (forward an backward)
%--------------------------------------------------------------------------
A{1} = [0 0; 0 0];
A{2} = [0 0; 0 0];
A{3} = [0 0; 0 0];
B    = {};
C    = sparse(2,0);
 
% get priors
%--------------------------------------------------------------------------
pE    = spm_dcm_neural_priors(A,B,C,options.model);
pE    = spm_L_priors(M.dipfit,pE);
pE    = spm_ssr_priors(pE);
[x,f] = spm_dcm_x_neural(pE,options.model);
 
% create forward model
%--------------------------------------------------------------------------
M.f   = f;
M.g   = 'spm_gx_erp';
M.x   = x;
M.n   = length(spm_vec(x));
M.pE  = pE;
M.m   = Ns;
M.l   = Nc;
M.Hz  = Hz;
M.Rft = 4;

% specify M.u - endogenous input (fluctuations) and intial states
%--------------------------------------------------------------------------
M.u   = sparse(Ns,1);
M.x   = spm_dcm_neural_x(pE,M);

% get expected CSD (observations)
%==========================================================================

% (log) connectivity parameters (forward connection only)
%--------------------------------------------------------------------------
pE.A{1}(2,1) = 2;
pE.S         = 1/8;

% (log) amplitude of fluctations and noise
%--------------------------------------------------------------------------
pE.a(1,:) = -2;
pE.b(1,:) = -4;
pE.c(1,:) = -4;

% expected cross spectral density
%--------------------------------------------------------------------------
[csd,Hz,mtf] = spm_csd_mtf(pE,M);

% place in data structure
%--------------------------------------------------------------------------
DCM.xY.y  = csd;
DCM.xY.dt = dt;
DCM.xY.Hz = Hz;


% Invert using a null model
%==========================================================================
DCM.options.model   = 'NULL';
DCM.options.spatial = 'LFP';
DCM.options.DATA    = 0;

DCM.M.dipfit.Nc  = Nc;
DCM.M.dipfit.Ns  = Ns;
DCM.M.U  = eye(Nc,Nc);

DCM.A = {ones(Nc,Nc)};
DCM.B = {};
DCM.C = sparse(Nc,0);

% estimate
%--------------------------------------------------------------------------
DCM   = spm_dcm_csd(DCM);


% show results in terms of transfer functions and Granger causality
%==========================================================================
spm_figure('GetWin','Figure 1'); clf

% transfer functions in the absence of measurement noise
%--------------------------------------------------------------------------
Ep    = DCM.Ep; 
Ep.b  = Ep.b - 32;              % and suppress non-specific and
Ep.c  = Ep.c - 32;              % specific channel noise

Gu           = spm_csd_mtf_gu(pE,Hz);
[psd Hz dtf] = spm_csd_mtf(Ep,DCM.M);

mtf   = mtf{1};
csd   = csd{1};
psd   = psd{1};


tew   = spm_dtf2gew(mtf,Gu);
ccf   = spm_csd2ccf(psd,Hz,dt);
qew   = spm_ccf2gew(ccf,Hz,dt,p);
ccf   = spm_csd2ccf(csd,Hz,dt);
gew   = spm_ccf2gew(ccf,Hz,dt,p);

spm_spectral_plot(Hz,tew, 'b', 'frequency','density')
spm_spectral_plot(Hz,qew, 'r',  'frequency','density')
spm_spectral_plot(Hz,gew, 'g',  'frequency','density')
legend('Granger causality (true)',...
       'Granger causality (source)',...
       'Granger causality (channel)')

subplot(2,2,3), a = axis;
subplot(2,2,1), axis(a);
subplot(2,2,2), axis(a);
subplot(2,2,4), axis(a);
