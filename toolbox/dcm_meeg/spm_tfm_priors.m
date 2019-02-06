function [E,V] = spm_tfm_priors(A,B,C)
% prior moments for a canonical microcircuit model (with plasticity)
% FORMAT [pE,pC] = spm_tfm_priors(A,B,C)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
%
% pE - prior expectation - f(x,u,P,M)
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - syaptic time constants
%    pE.S - intrinsic again
%    pE.G - intrinsic connection strengths
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A - extrinsic
%    pE.B - trial-dependent (driving)
%    pE.C - stimulus input
%    pE.D - delays
%    pE.N - trial-dependent (modulatory)
%
% plasticity parameters
%--------------------------------------------------------------------------
%    pE.E - voltage-dependent potentiation
%    pE.F - decay
%
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R - onset and dispersion
%
% pC - prior (co)variances
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_fx_cmc_tfm
%__________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%__________________________________________________________________________
% Copyright (C) 2011 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_tfm_priors.m 7527 2019-02-06 19:12:56Z karl $
 
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
n     = size(C,1);                                % number of sources
u     = size(C,2);                                % number of inputs
 
% parameters for neural-mass forward model
%==========================================================================
 
% restructure adjacency matrices
%--------------------------------------------------------------------------
D{1}  = A{1};                                     % forward  (i)
D{2}  = A{1};                                     % forward  (ii)
D{3}  = A{2};                                     % backward (i)
D{4}  = A{2};                                     % backward (ii)

% modulatory extrinsic connectivity
%--------------------------------------------------------------------------
if numel(A) < 3, A{3} = spm_zeros(A{1}); end

E.M   = 0*A{3};
V.M   = ~~A{3};
A     = D;
 
% extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
      A{i} = ~~A{i};
    E.A{i} = A{i}*32 - 32;
    V.A{i} = A{i}/8;
    Q      = Q | A{i};
end

% allow intrinsic delays
%--------------------------------------------------------------------------
Q     = Q | speye(n,n);                            
 
% driving connectivity - input-dependent scaling
%--------------------------------------------------------------------------
for i = 1:length(B)
      B{i} = ~~B{i};
    E.B{i} = 0*B{i};
    try
        V.B{i} = (B{i} & Q & ~V.M)/8;
    catch
        V.B{i} = B{i}/8;
    end
end

% plasticity parameters
%--------------------------------------------------------------------------
E.E  = sparse(n,4);   V.E  = sparse(n,4) + 1/16;  % potentiation                                             % potentiation
E.F  = sparse(n,4);   V.F  = sparse(n,4) + 1/16;  % decay


% modulatory connectivity - input-dependent scaling
%--------------------------------------------------------------------------
for i = 1:length(B)
    E.N{i} = 0*B{i};
    try
        V.N{i} = (B{i} & Q & V.M)/8;
    catch
        V.N{i} = spm_zeros(B{i});
    end
end

% exogenous connectivity - where inputs enter
%--------------------------------------------------------------------------
C    = ~~C;
E.C  = C*32 - 32;
V.C  = C/32;
 
% synaptic parameters
%--------------------------------------------------------------------------
E.T  = sparse(n,4);   V.T  = sparse(n,4) + 1/16;  % time constants
E.G  = sparse(n,4);   V.G  = sparse(n,4) + 1/16;  % intrinsic connectivity
E.D  = sparse(n,n);   V.D  = Q/32;                % delay
E.S  = 0;             V.S  = 1/16;                % intrinsic gain

 
% set stimulus parameters: onset and dispersion
%--------------------------------------------------------------------------
E.R    = sparse(u,2);  V.R   = ones(u,1)*[1/16 1/16];
warning('on','MATLAB:log:logOfZero');
 
return
 
 
 
% demo for log-normal pdf
%==========================================================================
x  = (1:64)/16;
for i = [2 16 32 128 256]
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
axis square
