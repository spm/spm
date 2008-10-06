function [varargout] = spm_sep_priors(A,B,C,dipfit)
% prior moments for a neural-mass model of erps
% FORMAT [pE,gE,pC,gC] = spm_sep_priors(A,B,C,dipfit)
% FORMAT [pE,pC]       = spm_sep_priors(A,B,C,dipfit)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
% dipfit       - forward model structure
%
% pE - prior expectation - f(x,u,P,M)
% gE - prior expectation - g(x,u,G,M)
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - synaptic time constants
%    pE.H - synaptic densities
%    pE.S - activation function parameters
%
% spatial parameters
%--------------------------------------------------------------------------
%    gE.Lpos - position                   - ECD
%    gE.L    - moment (orientation)       - ECD
%
% or gE.L    - coefficients of local modes - Imaging
%    gE.L    - gain of electrodes         - LFP
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A - extrinsic
%    pE.B - trial-dependent
%    pE.C - stimulus input
%    pE.G - intrinsic
%    pE.D - delays
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R - magnitude, onset and dispersion
%    pE.N - background fluctuations
%
% pC - prior covariances: cov(spm_vec(pE))
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_erp_fx
%
% This is like spm_erp_priors but:
%
%  (i) treats intrinsicn connection strneghts as free parapters
%  (ii) has fastr prior time constants
%  (iii) Allows ecah subpouaiton to express its own dipole
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_sep_priors.m 2310 2008-10-06 19:20:45Z karl $
 
% default: a single source model
%--------------------------------------------------------------------------
if nargin < 3
    A   = {0 0 0};
    B   = {};
    C   = 1;
end
 
% disable log zero warning
%--------------------------------------------------------------------------
warning('off','MATLAB:log:logOfZero');

n     = size(C,1);                                   % number of sources
u     = size(C,2);                                   % number of inputs
n1    = ones(n,1);
 
 
% parameters for electromagnetic forward model 
%==========================================================================
 
% source contribution to ECD of each state (1 7 9)
%--------------------------------------------------------------------------
try,   type = dipfit.type; catch, type = 'LFP'; end
for i = 1:3
    switch type
 
        case{'ECD'}
            %--------------------------------------------------------------
            G.Lpos = dipfit.Lpos;    U.Lpos =      0*ones(3,n); % position
            G.L{i} = sparse(3,n);    U.L{i} = exp(8)*ones(3,n); % moments
 
        case{'IMG'}
            %--------------------------------------------------------------
            m      = dipfit.Nm;
            G.L{i} = sparse(m,n);    U.L{i} = exp(8)*ones(m,n); % modes
 
        case{'LFP'}
            %--------------------------------------------------------------
            G.L{i} = ones(1,n);      U.L{i} = exp(8)*ones(1,n); % gains
            
    end
end
 
% parameters for neural-mass forward model
%==========================================================================
 
% set intrinsic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(n1) - 1;       V.T = n1/8;             % time constants
E.H   = log(n1);           V.H = n1/32;            % synaptic density
 
% set intrinsic [excitatory] connection strengths
%--------------------------------------------------------------------------
E.G   = n1*[0 0 0 0];      V.G = ones(n,4)/8;      % synaptic density
 
% set intrinsic [excitatory] time constants
%--------------------------------------------------------------------------
E.S   = [0 0];             V.S = [1 1]/8;          % dispersion & threshold
 
% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
      A{i} = ~~A{i};
    E.A{i} = A{i}*32 - 32;                         % forward
    V.A{i} = A{i}/2;                              % backward
    Q      = Q | A{i};                             % and lateral connections
end
 
for i = 1:length(B)
      B{i} = ~~B{i};
    E.B{i} = 0*B{i};                               % input-dependent scaling
    V.B{i} = B{i}/2;
    Q      = Q | B{i};
end
C      = ~~C;
E.C    = C*32 - 32;                                % where inputs enter
V.C    = C/32;
 
% set delay (enforcing symmetric delays)
%--------------------------------------------------------------------------
E.D        = sparse(n,n);
V.D        = Q/16;
 
% set stimulus parameters: onset and dispersion
%--------------------------------------------------------------------------
E.R        = ones(u,1)*[0 -2];  V.R  = ones(u,1)*[1/16 1/16];
 
warning('on','MATLAB:log:logOfZero');

 
% prior moments, if two arguments
%==========================================================================
if nargout == 4
    varargout{1} = E;
    varargout{2} = G;
    varargout{3} = diag(sparse(spm_vec(V)));
    varargout{4} = diag(sparse(spm_vec(U)));
    return
else
    switch type
 
        case{'ECD'}
            %--------------------------------------------------------------
            E.Lpos  = G.Lpos;
            E.L     = G.L   ;
            C       = diag(sparse(spm_vec(V,U)));
 
        case{'IMG','LFP'}
            %--------------------------------------------------------------
            E.L     = G.L;
            C       = diag(sparse(spm_vec(V,U)));
            
        otherwise
    end
 
end
 
if nargout == 2
    varargout{1} = E;
    varargout{2} = C;
    return
end
