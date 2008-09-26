function [varargout] = spm_erpsymm_priors(A,B,C,dipfit,ppC,pgC)
% prior moments for a neural-mass model of erps
% FORMAT [pE,gE,pC,gC] = spm_erp_priors(A,B,C,dipfit)
% FORMAT [pE,pC] = spm_erp_priors(A,B,C,dipfit)
% FORMAT [M]     = spm_erp_priors(A,B,C,dipfit)
% FORMAT           spm_erp_priors(A,B,C,dipfit)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
% dipfit       - prior forward model structure
%
% pE - prior expectation - f(x,u,P,M)
% gE - prior expectation - g(x,u,G,M)
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - syaptic time constants
%    pE.H - syaptic densities
%
% spatial parameters
%--------------------------------------------------------------------------
%    gE.Lpos - position                   - ECD
%    gE.L    - moment (orientation)       - ECD
%
% or gE.L    - coeficients of local modes - Imaging
%    gE.L    - gain of electrodes         - LFP
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A - extrinsic
%    pE.B - trial-dependent
%    pE.C - stimulus input
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
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_erpsymm_priors.m 2208 2008-09-26 18:57:39Z karl $

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

% paramters for electromagnetic forward model
%==========================================================================
try,   type = dipfit.type; catch, type = 'LFP'; end
switch type
    
    case{'ECD'}
    %----------------------------------------------------------------------
    G.Lpos = dipfit.Lpos;   U.Lpos =      0*ones(3,n);    % positions
    G.L    = sparse(3,n);   U.L    = exp(8)*ones(3,n);    % orientations
 
    case{'IMG'}
    %----------------------------------------------------------------------
    m      = dipfit.Nm;
    G.L    = sparse(m,n);   U.L    = exp(8)*ones(m,n);    % modes
    
    case{'LFP'}
    %----------------------------------------------------------------------
    G.L    = ones(1,n);     U.L    = exp(8)*ones(1,n);    % gains
    
end

% source contribution to ECD
%--------------------------------------------------------------------------
G.J   = sparse(1,[1 7 9],[0.2 0.2 0.6],1,9);
U.J   = G.J/128;

% parameters for neural-mass forward model
%==========================================================================

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(n1);        V.T = n1/8;                % time constants
E.H   = log(n1);        V.H = n1/8;                % synaptic density

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.S   = [0 0];          V.S = [1 1]/8;             % dispersion & threshold


% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
      A{i} = ~~A{i};
    E.A{i} = A{i}*32 - 32;                         % forward
    V.A{i} = A{i}/16;                              % backward
    Q      = Q | A{i};                             % and lateral connections
end
 
for i = 1:length(B)
      B{i} = ~~B{i};
    E.B{i} = 0*B{i};                               % input-dependent scaling
    V.B{i} = B{i}/16;
    Q      = Q | B{i};
end
C      = ~~C;
E.C    = C*32 - 32;                                % where inputs enter
V.C    = C/32;

% set delay (enforcing symmetric delays)
%--------------------------------------------------------------------------
E.D        = sparse(n,n);
V.D        = Q/8;

% set stimulus parameters: onset and dispersion
%--------------------------------------------------------------------------
E.R        = sparse(u,2);  V.R   = ones(u,1)*[1/16 1/16];

warning('on','MATLAB:log:logOfZero');


% prior momments if two arguments
%--------------------------------------------------------------------------
if nargout == 4
    varargout{1} = E;
    varargout{2} = G;
    
    % replace covariance matrices by carfully constructed input
    %----------------------------------------------------------------------
    varargout{3} = ppC;
    varargout{4} = pgC;
    return
else
    switch type

        case{'ECD'}
            %--------------------------------------------------------------
            E.Lpos  = G.Lpos;
            E.L     = G.L;
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
