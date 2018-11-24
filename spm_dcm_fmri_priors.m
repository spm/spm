function [pE,C,x,pC] = spm_dcm_fmri_priors(A,B,C,D,options)
% Returns the priors for a two-state DCM for fMRI.
% FORMAT:[pE,pC,x,vC] = spm_dcm_fmri_priors(A,B,C,D,options)
%
%   options.two_state:  (0 or 1) one or two states per region
%   options.stochastic: (0 or 1) exogenous or endogenous fluctuations
%   options.precision:           log precision on connection rates
%
% INPUT:
%    A,B,C,D - constraints on connections (1 - present, 0 - absent)
%
% OUTPUT:
%    pE     - prior expectations (connections and hemodynamic)
%    pC     - prior covariances  (connections and hemodynamic)
%    x      - prior (initial) states
%    vC     - prior variances    (in struct form)
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
% $Id: spm_dcm_fmri_priors.m 7497 2018-11-24 17:00:25Z karl $

% number of regions
%--------------------------------------------------------------------------
n = length(A);

% check options and D (for nonlinear coupling)
%--------------------------------------------------------------------------
try, options.stochastic; catch, options.stochastic = 0; end
try, options.induced;    catch, options.induced    = 0; end
try, options.two_state;  catch, options.two_state  = 0; end
try, options.backwards;  catch, options.backwards  = 0; end
try, D;                  catch, D = zeros(n,n,0);       end


% connectivity priors and intitial states
%==========================================================================
if options.two_state
    
    % (6) initial states
    %----------------------------------------------------------------------
    x     = sparse(n,6);
    A     = logical(A - diag(diag(A)));
    
    % precision of log-connections (two-state)
    %---------------------------------------------------------------------
    try, pA = exp(options.precision); catch,  pA = 16;  end
    
    % prior expectations and variances
    %----------------------------------------------------------------------
    pE.A  =  (A + eye(n,n))*32 - 32;
    pE.B  =  B*0;
    pE.C  =  C*0;
    pE.D  =  D*0;
    
    % prior covariances
    %----------------------------------------------------------------------
    for i = 1:size(A,3)
        pC.A(:,:,i) = A(:,:,i)/pA + eye(n,n)/pA;
    end
    pC.B  =  B/4;
    pC.C  =  C*4;
    pC.D  =  D/4;
    
    % excitatory proportion
    %----------------------------------------------------------------------
    if options.backwards
        pE.A(:,:,2) = A*0;
        pC.A(:,:,2) = A/pA;
    end

else
    
    % one hidden state per node
    %======================================================================
    
    % (6 - 1) initial states
    %----------------------------------------------------------------------
    x     = sparse(n,5);
    
    % precision of connections (one-state)
    %---------------------------------------------------------------------
    try, pA = exp(options.precision); catch,  pA = 64;  end
    try, dA = options.decay;          catch,  dA = 1;   end
    
    % prior expectations
    %----------------------------------------------------------------------
    if isvector(A)
        A     = logical(A);
        pE.A  = (A(:) - 1)*dA;
    else
        A     = logical(A - diag(diag(A)));
        pE.A  = A/128;
    end
    pE.B  = B*0;
    pE.C  = C*0;
    pE.D  = D*0;
    
    % prior covariances
    %----------------------------------------------------------------------
    if isvector(A)
        pC.A  = A(:);
    else
        for i = 1:size(A,3)
            pC.A(:,:,i) = A(:,:,i)/pA + eye(n,n)/pA;
        end
    end
    pC.B  = B;
    pC.C  = C;
    pC.D  = D;
    
end

% and add hemodynamic priors
%==========================================================================
pE.transit = sparse(n,1);  pC.transit = sparse(n,1) + 1/256;
pE.decay   = sparse(1,1);  pC.decay   = sparse(1,1) + 1/256;
pE.epsilon = sparse(1,1);  pC.epsilon = sparse(1,1) + 1/256;


% add prior on spectral density of fluctuations (amplitude and exponent)
%--------------------------------------------------------------------------
if options.induced
    pE.a  = sparse(2,1);   pC.a = sparse(2,1) + 1/64; % neuronal fluctuations
    pE.b  = sparse(2,1);   pC.b = sparse(2,1) + 1/64; % channel noise global
    pE.c  = sparse(1,n);   pC.c = sparse(1,n) + 1/64; % channel noise specific
end

% prior covariance matrix
%--------------------------------------------------------------------------
C  = diag(spm_vec(pC));

return
