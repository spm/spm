function spm_delays_demo
% Demo routine for induced responses
%==========================================================================
%
% This routine illustrates the Taylor approximation to delay differential
% equation solvers using two (extrinsically connected) neural masses. In
% this simulation, using a canonical microcircuit model, exogenous inputs
% are applied to two sources with a unidirectional (forward) connection.
% The responses of those regions are summarised in terms of their
% first-order Volterra kernels, under different conduction delays from the
% source to the target. The effect of these delays can then be seen as a
% translation of the forward curve and (or impulse response of the target 
% to perturbations of the source.
% 
% See also:
%  spm_dcm_delay.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_delays_demo.m 6900 2016-10-08 13:16:46Z karl $
 

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
pE.G(:,3)    = 0;

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
k      = (1:8)/8;
M.pst  = (1:128)/1000;
for j  = 1:length(k)
    
    % keep total power of fluctuations constant
    %----------------------------------------------------------------------
    P        = pE;
    P.D(2,1) = log(k(j));

    % create forward model and solve for steady state
    %----------------------------------------------------------------------
    M.x      = spm_dcm_neural_x(P,M);
    
    % first-order Volterra kernels
    %======================================================================
    [S,K,s,w,t]  = spm_dcm_mtf(P,M);
    
    spm_spectral_plot(t*1000,K,'r','peristimulus time (ms)','density',1)
    
end

subplot(2,2,1); title('response of source','FontSize',16)
subplot(2,2,3); title('forward influence ','FontSize',16); a = axis;
subplot(2,2,2); title('backward influence','FontSize',16); axis(a) 
subplot(2,2,4); title('response of target','FontSize',16)

