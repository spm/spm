function [pE,pC] = spm_mfm_priors(A,B,C,L,J)
% prior moments for a neural mass model of erps
% FORMAT [pE,pC] = spm_mfm_priors(A,B,C,L,J)
%
% A{3},B,C    - binary constraints on extrinsic connectivity
% L           - Lead-field parameters
%
% pE - prior expectation
% pC - prior covariance
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - synaptic time constants
%    pE.H - synaptic densities
%
% spatial parameters
%--------------------------------------------------------------------------
%    pE.L - moment (orientation), mode or gain parameters; plus
%    pE.Lpos - position for ECD models
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A - extrinsic - coupling
%    pE.B - extrinsic - input-dependent
%    pE.C - extrinsic - stimulus input
%    pE.G - intrinsic
%
%--------------------------------------------------------------------------
%
% pC - prior covariances: cov(spm_vec(pE))
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_lfp_fx or M.fP
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_mfm_priors.m 1174 2008-02-27 20:22:30Z karl $

% defaults
%--------------------------------------------------------------------------
if nargin <  3                                      % a single source model
    A   = {0 0 0};
    B   = {};
    C   = 1;
end
n   = size(C,1);                                    % number of sources

% disable log zero warning
%--------------------------------------------------------------------------
warning off


% paramters for electromagnetic forward model
%==========================================================================

% source-specific contribution to LFP
%--------------------------------------------------------------------------
try, J; catch, J  = sparse([1 7 9],1,[0.2 0.2 0.6],13,1); end

% source-specific contribution to LFP
%--------------------------------------------------------------------------
E.J   = J;                 V.J = J/128;            % contributing states

% lead field parameters
%--------------------------------------------------------------------------
try,   type = L.type; catch, type = 'LFP'; end
switch type
    
    case{'ECD (EEG)','ECD (MEG)'}
    %----------------------------------------------------------------------
    E.Lpos = L.L.pos;       V.Lpos =   0*ones(3,n);    % positions
    E.L    = sparse(3,n);   V.L    = 256*ones(3,n);    % orientations

    case{'Imaging'}
    %----------------------------------------------------------------------
    m      = L.Nm;
    E.L    = sparse(m,n);   V.L    =  16*ones(m,n);    % modes
    
    case{'LFP'}
    %----------------------------------------------------------------------
    E.L    = ones(1,n);     V.L    = 256*ones(1,n);    % gains
    
end


% parameters for neural-mass forward model
%==========================================================================

% sigmoid parameters
%--------------------------------------------------------------------------
E.R   = [0 0];             V.R = [1 1]/8;

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(ones(n,2));    V.T = ones(n,2)/8;      % time constants
E.G   = log(ones(n,2));    V.G = ones(n,2)/8;      % intrinsic connections

% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
    E.A{i} = log(~~A{i} + eps);                    % forward
    V.A{i} = A{i}/2;                               % backward
    Q      = Q | A{i};                             % and lateral connections
end

for i = 1:length(B)
    E.B{i} = 0*B{i};                               % input-dependent scaling
    V.B{i} = ~~B{i}/2;
    Q      = Q | B{i};
end
E.C    = log(~~C + eps);                           % where inputs enter
V.C    = C/2;

% vectorize
%--------------------------------------------------------------------------
pE     = E;
pV     = spm_vec(V);
pC     = diag(sparse(pV));
warning on

return



