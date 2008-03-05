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
%    pE.T    - syaptic time constants
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
% $Id: spm_nmm_priors.m 1190 2008-03-05 17:20:23Z karl $
 
% default: a single source model
%--------------------------------------------------------------------------
if nargin < 3
    A   = {0 0 0};
    B   = {};
    C   = 1;
end
 
% disable log zero warning
%--------------------------------------------------------------------------
warning off
n     = size(C,1);                                   % number of sources
u     = size(C,2);                                   % number of inputs
 
% paramters for electromagnetic forward model
%==========================================================================
try,   type = dipfit.type; catch, type = 'LFP'; end
switch type
    
    case{'ECD (EEG)','ECD (MEG)'}
    %----------------------------------------------------------------------
    G.Lpos = dipfit.L.pos;  U.Lpos =      0*ones(3,n);    % positions
    G.L    = sparse(3,n);   U.L    = exp(8)*ones(3,n);    % orientations
 
    case{'Imaging'}
    %----------------------------------------------------------------------
    m      = dipfit.Nm;
    G.L    = sparse(m,n);   U.L    = exp(8)*ones(m,n);    % modes
    
    case{'LFP'}
    %----------------------------------------------------------------------
    G.L    = ones(1,n);     U.L    = exp(8)*ones(1,n);    % gains
    
end
 
% contribution of vec(x) to ECD (first 3 populations)
%--------------------------------------------------------------------------
G.J   = [0.2 0.2 0.6];
U.J   = G.J/128;
 
% parameters for neural-mass forward model
%==========================================================================
 
% population variance
%--------------------------------------------------------------------------
E.S   = 0;                 V.S = 1/8;

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(ones(n,1));    V.T = ones(n,1)/8;      % excitatory constants
E.G   = log(ones(n,1));    V.G = ones(n,1)/8;      % intrinsic connections
 
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
C      = ~~C;
E.C    = C*32 - 32;                                % where inputs enter
V.C    = C;

 
% set stimulus parameters: onset and dispersion
%--------------------------------------------------------------------------
E.R        = sparse(u,2);  V.R   = ones(u,1)*[1/16 1/16];
warning on
 
% prior moments if two arguments
%==========================================================================
if nargout == 4
    varargout{1} = E;
    varargout{2} = G;
    varargout{3} = diag(sparse(spm_vec(V)));
    varargout{4} = diag(sparse(spm_vec(U)));
    return
else
    switch type
 
        case{'ECD (EEG)','ECD (MEG)'}
            %--------------------------------------------------------------
            E.Lpos  = G.Lpos;
            E.L     = G.L;
            E.J     = G.J;
            C       = diag(sparse(spm_vec(V,U)));
 
        case{'Imaging','LFP'}
            %--------------------------------------------------------------
            E.L     = G.L;
            E.J     = G.J;
            C       = diag(sparse(spm_vec(V,U)));
            
        otherwise
    end
 
end
 
if nargout == 2
    varargout{1} = E;
    varargout{2} = C;
    return
end
