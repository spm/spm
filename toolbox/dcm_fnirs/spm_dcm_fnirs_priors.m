function [pE,pC,x] = spm_dcm_fnirs_priors(DCM)
% Returns the priors for a one-state DCM for fNIRS.
% FORMAT:[pE,pC,x] = spm_dcm_fnirs_priors(DCM)
%
% INPUT:
%    DCM.a,DCM.b,DCM.c,DCM.c - constraints on connections (1 - present, 0 - absent)
%    DCM.n - number of sources of interest 
%    DCM.Y.nch - number of channels of interest 
%    DCM.options.two_state:  (0 or 1) one or two states per region
%    DCM.options.stochastic: (0 or 1) exogenous or endogenous fluctuations
%    DCM.options.precision:           log precision on connection rates
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
%
% 3. Tak S, Kempny AM, Friston, KJ, Leff, AP, Penny WD. Dynamic causal
%    modelling for functional near-infrared spectroscopy. 
%    Neuroimage 111: 338-349, 2015. 
%
% This script is based on spm_dcm_fmri_priors.m written by Karl Friston.
% 
% In this script, optics priors are added, prior covariance of A is changed, 
% prior for extended Balloon model (viscoelastic time constant) is added. 
%__________________________________________________________________________
% Copyright (C) 2015 Wellcome Trust Centre for Neuroimaging

% Will Penny & Sungho Tak
% $Id: spm_dcm_fnirs_priors.m 6942 2016-11-21 13:17:44Z guillaume $

% Unpack
%--------------------------------------------------------------------------
A = DCM.a; 
B = DCM.b; 
C = DCM.c; 
D = DCM.d; 
n = DCM.n; % number of sources of interest 
nch = DCM.Y.P.nch; % number of channels of interest 
options = DCM.options; 

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
    x     = sparse(n,7);
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
    x     = sparse(n,6);
    
    % precision of connections (one-state) 
    %---------------------------------------------------------------------
    try, pA = exp(options.precision); catch,  pA = 16;  end
    try, dA = options.decay;          catch,  dA = 1;   end
    
    % prior expectations
    %----------------------------------------------------------------------
    if isvector(A)
        A     = logical(A);
        pE.A  = (A(:) - 1)*dA;
    else
        A     = logical(A - diag(diag(A)));
        pE.A  =  A/128;
    end
    pE.B  =  B*0;
    pE.C  =  C*0;
    pE.D  =  D*0;
    
    % prior covariances
    %----------------------------------------------------------------------
    if isvector(A)
        pC.A  = A(:);
    else
        for i = 1:size(A,3)
            pC.A(:,:,i) = A(:,:,i)/pA + eye(n,n)/pA;
        end
    end
    pC.B  =  B;
    pC.C  =  C;
    pC.D  =  D;
    
end

% add hemodynamic priors
%--------------------------------------------------------------------------
pE.transit = sparse(n,1);  pC.transit = sparse(n,1) + exp(-6);
pE.decay   = sparse(n,1);  pC.decay   = sparse(n,1) + exp(-6); 
pE.afback = sparse(n,1); pC.afback = sparse(n,1) + exp(-6); % transit time 
pE.tv = sparse(n,1); pC.tv = sparse(n,1) + 1; % viscoelastic time 

% add optics priors
%--------------------------------------------------------------------------
if DCM.options.pialv % correction of pial vein oxygenation
    pE.pv = sparse(nch.*2, 1); pC.pv = sparse(nch.*2, 1) + exp(-3); 
else
    pE.pv = sparse(nch, 1); pC.pv = sparse(nch, 1) + exp(-3); 
end

if options.rs ~= 0 % spatially distributed sources 
    pE.rs = sparse(n,1); pC.rs = sparse(n,1) + exp(-6); 
end

% add prior on spectral density of fluctuations (amplitude and exponent)
%--------------------------------------------------------------------------
if options.induced
    pE.a  = sparse(2,n);   pC.a = sparse(2,n) + 1/64; % neuronal fluctuations
    pE.b  = sparse(2,1);   pC.b = sparse(2,1) + 1/64; % channel noise global
    pE.c  = sparse(2,n);   pC.c = sparse(2,n) + 1/64; % channel noise specific
end

% prior covariance matrix
%--------------------------------------------------------------------------
pC  = diag(spm_vec(pC));

return
