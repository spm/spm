function [varargout] = spm_nmm_priors(A,B,C,dipfit)
% prior moments for a neural-mass model of ERPs
% FORMAT [pE,gE,pC,gC] = spm_nmm_priors(A,B,C,dipfit)
% FORMAT [pE,pC]       = spm_nmm_priors(A,B,C,dipfit)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
% dipfit       - prior forward model structure
%
% pE - prior expectation - f(x,u,P,M)
% gE - prior expectation - g(x,u,G,M)
%
% spatial parameters
%--------------------------------------------------------------------------
%    gE.Lpos - position                   - ECD
%    gE.L    - orientation
%
% or gE.L    - coefficients of local modes - Imaging
%    gE.L    - gain of electrodes         - LFP
%
% source contribution to ECD
%--------------------------------------------------------------------------
%    gE.J    - contribution
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T    - synaptic time constants
%    pE.G    - intrinsic connectivity
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A    - extrinsic
%    pE.B    - trial-dependent
%    pE.C    - stimulus input
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R    - onset and dispersion
%
% pC - prior covariances: cov(spm_vec(pE))
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_erp_fx
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nmm_priors.m 1277 2008-03-28 18:36:49Z karl $
 
 
% disable log zero warning
%--------------------------------------------------------------------------
warning('off','MATLAB:log:logOfZero');
n     = size(C,1);                                   % number of sources
u     = size(C,2);                                   % number of inputs
 
 
% parameters for electromagnetic forward model
%==========================================================================
try
    [G,U] = spm_L_priors(dipfit);
catch
    [G,U] = spm_L_priors(n);
end

% contribution of states to ECD (first 3 populations)
%--------------------------------------------------------------------------
G.J   = [0.1 0.1 1];
U     = spm_cat(diag({U, diag([1/128 1/512 0])}));
 


% parameters for neural-mass forward model
%==========================================================================
 
% population variance
%--------------------------------------------------------------------------
E.S   = 0;                 V.S = 1/64;
 
% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(ones(n,1));    V.T = ones(n,1)/64;      % excitatory constants
E.G   = log(ones(n,1));    V.G = ones(n,1)/64;      % intrinsic connections
 
% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
    A{i} = ~~A{i};
    E.A{i} = A{i}*32 - 32;                         % forward
    V.A{i} = A{i};                                 % backward
    Q      = Q | A{i};                             % and lateral connections
end
 
for i = 1:length(B)
    B{i} = ~~B{i};
    E.B{i} = 0*B{i};                               % input-dependent scaling
    V.B{i} = B{i};
    Q      = Q | B{i};
end
C     = ~~C;
E.C   = C*33 - 32;                                 % where inputs enter
V.C   = C;
 
 
% set stimulus parameters: onset and dispersion
%--------------------------------------------------------------------------
E.R   = sparse(u,2);  V.R   = ones(u,1)*[1/16 1/16];
warning('on','MATLAB:log:logOfZero');
 
 
% covariance of neural parameters
%==========================================================================
V     = diag(sparse(spm_vec(V)));

% format output arguments
%--------------------------------------------------------------------------
if nargout == 4
    
    varargout{1} = E;
    varargout{2} = G;
    varargout{3} = V;
    varargout{4} = U;
 
elseif nargout == 2
 
    % transcribe spatial parameters
    %----------------------------------------------------------------------
    E.Lpos  = G.Lpos;
    E.L     = G.L;
    E.J     = G.J;
 
    varargout{1} = E;
    varargout{2} = spm_cat(diag({V,U}));
    
end
