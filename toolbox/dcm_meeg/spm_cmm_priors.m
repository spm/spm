function [pE,pC] = spm_cmm_priors(A,B,C)
% prior moments for a canonical neural-mass model of ERPs
% FORMAT [pE,pC] = spm_cmm_priors(A,B,C)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
%
% pE - prior expectation - f(x,u,P,M)
%
% population variance
%--------------------------------------------------------------------------
%     E.S    - variance
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
%    pE.D    - delays
%    pE.U    - exogenous background activity
%
% pC - prior (co)variances
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
% $Id: spm_cmm_priors.m 4814 2012-07-30 19:56:05Z karl $
 
 
% disable log zero warning
%--------------------------------------------------------------------------
warning('off','MATLAB:log:logOfZero');
n     = size(C,1);                                % number of sources
u     = size(C,2);                                % number of inputs
p     = 4;                                        % number of populations
 
% parameters for neural-mass forward model
%==========================================================================
 
% population variance
%--------------------------------------------------------------------------
pE.S  = 0;                 pC.S = 1/16;
 
% intrinic [excitatory] time constants (H mediates the effects of B)
%--------------------------------------------------------------------------
pE.T  = zeros(n,1);        pC.T =  ones(n,1)/16;
pE.H  = zeros(n,1);        pC.H = zeros(n,1)/16;


% extrinsic connectivity (n x n)
%==========================================================================

% restructure adjacency matrices
%--------------------------------------------------------------------------
D{1}  = A{1};                                     % forward  (i)
D{2}  = A{1};                                     % forward  (ii)
D{3}  = A{2};                                     % backward (i)
D{4}  = A{2};                                     % backward (ii)
A     = D;
Q     = sparse(n,n);
for i = 1:length(A)
    A{i}    = ~~A{i};
    pE.A{i} = A{i}*32 - 32;                       % forward
    pC.A{i} = A{i}/8;                             % backward
    Q       = Q | A{i};                           % and lateral connections
end
 
% input-dependent scaling
%--------------------------------------------------------------------------
for i = 1:length(B)
    B{i} = ~~B{i};
    pE.B{i} = 0*B{i};
    pC.B{i} = B{i}/8;
    Q       = Q | B{i};
end

% exogenous inputs
%--------------------------------------------------------------------------
C     = ~~C;
pE.C  = C*32 - 32;
pC.C  = C/32;


% intrinsic connectivity (p x p x n)
%==========================================================================
pE.G  = repmat(zeros(p,p),[1 1 n]);
pC.G  = repmat(  eye(p,p),[1 1 n]);


% Exogenous inputs: onset and dispersion
%--------------------------------------------------------------------------
pE.R  = zeros(u,2);     pC.R  = ones(u,1)*[1/16 1/16];

% Delays (intrinsic and extrinsic)
%--------------------------------------------------------------------------
pE.D  = [0 0];          pC.D  = [1 1]/64;

% Capacitance
%--------------------------------------------------------------------------
pE.CV = 0;              pC.CV = 1/16;


 


warning('on','MATLAB:log:logOfZero');
