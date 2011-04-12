function [E, V ] = spm_nfm_priors(A,B,C)
% prior moments for a neural mass model of ERPs
% FORMAT [pE,pC] = spm_nfm_priors(A,B,C)
%
% A{3},B{m},C    - binary constraints on extrinsic connectivity
%
% pE - prior expectation
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - synaptic time constants
%    pE.H - synaptic densities
%    pE.R - activation function parameters
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A - extrinsic - coupling
%    pE.B - extrinsic - trial-dependent
%    pE.C - extrinsic - stimulus input
%    pE.G - intrinsic
%    pE.D - extrinsic delays
%    pE.I - intrinsic delays
%
% spatial parameters
%--------------------------------------------------------------------------
%   pE.eps - inverse velocity
%   pE.ext - dispersion 
%   pE.A31 ]
%   pE.A12 ]  coupling parameters - single source 
%   pE.A31 ]
%
%--------------------------------------------------------------------------
% pC - prior covariances: cov(spm_vec(pE))
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_fx_erp_nfs2
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_nfm_priors.m 4305 2011-04-12 18:15:32Z karl $
 
% defaults
%--------------------------------------------------------------------------
if nargin < 3                                       % a single source model
    A   = {0 0 0};
    B   = {0};
    C   = 1;
end
n   = size(C,1);                                    % number of sources
 
% disable log zero warning
%--------------------------------------------------------------------------
warning('off','MATLAB:log:logOfZero');
 
% parameters for neural-field forward model
%==========================================================================
 
% sigmoid parameters
%--------------------------------------------------------------------------
E.R   = 0;                 V.R = 1/64;
 
% set intrinsic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(ones(n,1));    V.T = ones(n,1)/16;     % time constants
E.H   = log(ones(n,1));    V.H = ones(n,1)/16;     % synaptic density

% set intrinsic connections
%--------------------------------------------------------------------------
E.G   = log(ones(n,4));    V.G = ones(n,4)/16;     % intrinsic connections
 
% set spatial parameters
%--------------------------------------------------------------------------
E.vel = 0;                  V.vel = 1/16;
E.ext = 0;                  V.ext = 1/16;

% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
      A{i} = ~~A{i};
    E.A{i} = A{i}*32 - 32;                        % forward
    V.A{i} = A{i}/16;                             % backward
    Q      = Q | A{i};                            % and lateral connections
end
 
for i = 1:length(B)
      B{i} = ~~B{i};
    E.B{i} = 0*B{i};                              % input-dependent scaling
    V.B{i} = B{i}/8;
    Q      = Q | B{i};
end
C      = ~~C;
E.C    = C*32 - 32;                               % where inputs enter
V.C    = C/32;
 

warning('on','MATLAB:log:logOfZero');

