function [pE,pC,x] = spm_dcm_fmri_priors(A,B,C,D,options)
% Returns the priors for a two-state DCM for fMRI.
% FORMAT:[pE,pC,x] = spm_dcm_fmri_priors(A,B,C,D,options)
%
%   options.two_state:  (0 or 1) one or two states per region
%   options.stochastic: (0 or 1) exogenous or endogenous fluctuations
%
% INPUT:
%    A,B,C,D - constraints on connections (1 - present, 0 - absent)
%
% OUTPUT:
%    pE     - prior expectations (connections and hemodynamic)
%    pC     - prior covariances  (connections and hemodynamic)
%    x      - prior (initial) states
%__________________________________________________________________________
%
% References for state equations:
% 1. Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
%    fMRI: a two-state model.
%    Neuroimage. 2008 Jan 1;39(1):269-78.
%
% 2. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_fmri_priors.m 5601 2013-08-11 21:19:53Z karl $

% number of regions
%--------------------------------------------------------------------------
n = length(A);

% check options and D (for nonlinear coupling)
%--------------------------------------------------------------------------
try, options.two_state;  catch, options.two_state  = 0; end
try, D;                  catch, D = zeros(n,n,0);       end


% prior (initial) states and shrinkage priors on A for endogenous DCMs
%--------------------------------------------------------------------------
if options.two_state,  x = sparse(n,6); else, x = sparse(n,5); end


% connectivity priors and intitial states
%==========================================================================
if options.two_state
    
    % (6) initial states
    %----------------------------------------------------------------------
    x     = sparse(n,6);
    
    % enforce optimisation of intrinsic (I to E) connections
    %----------------------------------------------------------------------
    A     = (A + eye(n,n)) > 0;

    % prior expectations and variances
    %----------------------------------------------------------------------
    pE.A  =  A*32 - 32;
    pE.B  =  B*0;
    pE.C  =  C*0;
    pE.D  =  D*0;

    % prior covariances
    %----------------------------------------------------------------------
    pC.A  =  A/4;
    pC.B  =  B/4;
    pC.C  =  C*4;
    pC.D  =  D/4;

else
    
    % (6 - 1) initial states
    %----------------------------------------------------------------------
    x     = sparse(n,5);

    % enforce self-inhibition
    %----------------------------------------------------------------------
    A     =  A > 0;
    A     =  A - diag(diag(A));

    % prior expectations
    %----------------------------------------------------------------------
    pE.A  =  A/(64*n);
    pE.B  =  B*0;
    pE.C  =  C*0;
    pE.D  =  D*0;
    
    % prior covariances
    %----------------------------------------------------------------------
    pC.A  =  A*8/n + eye(n,n)/(64*n);
    pC.B  =  B;
    pC.C  =  C;
    pC.D  =  D;

end

% and add hemodynamic priors
%==========================================================================
pE.transit = sparse(n,1);  pC.transit = sparse(n,1) + exp(-6);
pE.decay   = sparse(n,1);  pC.decay   = sparse(n,1) + exp(-6);
pE.epsilon = sparse(1,1);  pC.epsilon = sparse(1,1) + exp(-6);


% add prior on spectral density of fluctuations (amplitude and exponent)
%--------------------------------------------------------------------------
if strcmp(options.analysis,'CSD')
    
    pE.a  = sparse(2,n);   pC.a = sparse(2,n) + 1/16; % neuronal fluctuations
    pE.b  = sparse(2,1);   pC.b = sparse(2,1) + 1/16; % channel noise global
    pE.c  = sparse(2,n);   pC.c = sparse(2,n) + 1/16; % channel noise specific
    
end

% prior covariance matrix
%--------------------------------------------------------------------------
pC         = diag(spm_vec(pC));

return
