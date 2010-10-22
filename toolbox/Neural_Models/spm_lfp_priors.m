function [E,V] = spm_lfp_priors(A,B,C)
% prior moments for a neural mass model of ERPs
% FORMAT [pE,pC] = spm_lfp_priors(A,B,C)
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
%--------------------------------------------------------------------------
%
% pC - prior (co)variances
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_lfp_fx
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_lfp_priors.m 4096 2010-10-22 19:40:34Z karl $
 
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
 
% parameters for neural-mass forward model
%==========================================================================
 
% sigmoid parameters
%--------------------------------------------------------------------------
E.R   = [0 0];             V.R = [1 1]/8;
 
% set intrinsic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(ones(n,2));    V.T = ones(n,2)/8;      % time constants
E.H   = log(ones(n,1));    V.H = ones(n,1)/16;     % synaptic density

% set intrinsic connections
%--------------------------------------------------------------------------
E.G   = log(ones(n,5));    V.G = ones(n,5)/16;     % intrinsic connections
 
 
% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
      A{i} = ~~A{i};
    E.A{i} = log(A{i} + eps);                      % forward
    V.A{i} = A{i}/2;                               % backward
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
 
% set delay
%--------------------------------------------------------------------------
E.D    = sparse(n,n);     V.D = Q/16;              % extrinsic delays
E.I    = 0;               V.I = 1/32;              % intrinsic delays
warning('on','MATLAB:log:logOfZero');


return
 
% demo for log-normal pdf
%--------------------------------------------------------------------------
x    = [1:64]/16;
for i = [2 16]
    v = 1/i;
    p = 1./x.*exp(-log(x).^2/(2*v))/sqrt(2*pi*v);
    plot(x,p)
    text(x(16),p(16),sprintf('variance = 1/%i',1/v))
    hold on
end
xlabel('scaling')
ylabel('density')
grid on
hold off
