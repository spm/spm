function spm_delays_demo
% Demo routine for induced responses
%==========================================================================
%
% This routine illustrates the Taylor approxmiation to deay differential
% equatiton solvers using two (exttrinsically connected) neural masses
% 
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m, 
%  spm_csd2coh.m, spm_ccf2gew, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and 
%  spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_delays_demo.m 5883 2014-02-18 10:32:23Z karl $
 

% Notes: analysis of delay operator
%==========================================================================
spm_figure('GetWin','delays'); clf

% Model specification
%==========================================================================
rng('default')
 
% number of regions
%--------------------------------------------------------------------------
Nc    = 2;                                       % number of channels
Ns    = 2;                                       % number of sources
Hz    = 1:64;                                    % frequency
options.spatial  = 'LFP';
options.model    = 'CMC';
M.dipfit.model = options.model;
M.dipfit.type  = options.spatial;
M.dipfit.Nc    = Nc;
M.dipfit.Ns    = Ns;

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

% (log) connectivity parameters
%--------------------------------------------------------------------------
pE.A{1}(2,1) = 0;

% orders and model
%==========================================================================
nx    = length(spm_vec(x));
 
% create forward model
%--------------------------------------------------------------------------
M.f   = f;
M.g   = 'spm_gx_erp';
M.x   = x;
M.n   = nx;
M.pE  = pE;
M.m   = Ns;
M.l   = Nc;
M.Hz  = Hz;
M.u   = sparse(Ns,1);
M.x   = spm_dcm_neural_x(pE,M);


% delays
%--------------------------------------------------------------------------
D      = 12;
k      = log((1:8)/8);
M.pst  = (1:128)/1000;
M.pF.D = [1 D];
for j  = 1:length(k)
    
    % keep total power of fluctuations constant
    %----------------------------------------------------------------------
    P        = pE;
    P.D(2,1) = k(j);

    % create forward model and solve for steady state
    %----------------------------------------------------------------------
    M.x      = spm_dcm_neural_x(P,M);
    
    % first-order Volterra kernels
    %======================================================================
    [S,K,s,w,t]  = spm_dcm_mtf(P,M);
    
    spm_spectral_plot(t*1000,K,'r','frequency','density',1)
    
end

