function [varargout] = spm_erp_priors(A,B,C,L,U)
% prior moments for a neural mass model of erps
% FORMAT [pE,pC] = spm_erp_priors(A,B,C,L,U)
% FORMAT [M]     = spm_erp_priors(A,B,C,L,U)
% FORMAT           spm_erp_priors(A,B,C,L,U)
%
% A{3},B{m},C    - binary constraints on extrinsic connectivity
% L              - lead field matrix         [default = 1 s]
% U              - trial duration (sec)      [default = 1 s]
%
% pE - prior expectation
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - syaptic time constants
%    pE.H - syaptic densities
%
% spatial parameters
%--------------------------------------------------------------------------
%    pE.Lpos - position
%    pE.Lmon - moment (orientation):
%
% connectivity parameters
%--------------------------------------------------------------------------
%    pE.A - extrinsic
%    pE.B - trial-dependent
%    pE.C - stimulus input
%    pE.D - delays
%
% stimulus parameters
%--------------------------------------------------------------------------
%    pE.R - onset and dispersion (Gamma paramters)
%    pE.N - background fluctuation (DCT parameters)
%    pE.U - trial dutation (secs)
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
% %W% Karl Friston %E%

% defaults
%--------------------------------------------------------------------------
if nargin <  3                                 % a single source model
    A   = {0 0 0};
    B   = {};
    C   = 1;
end
if nargin <  4, L = speye(length(C)); end
if nargin <  5, U = 1;                end

warning off

% disable log zero warning
%--------------------------------------------------------------------------
n     = size(C,1);                                 % number of sources
m     = 8;                                         % and noise components
n1    = ones(n,1);

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(n1);        V.T = n1/16;               % time constants
E.H   = log(n1);        V.H = n1/16;               % synaptic density

% set observer parameters
%--------------------------------------------------------------------------
E.Lpos = L.pos;         V.Lpos = L.Vpos;           % dipole positions
E.Lmom = L.mom;         V.Lmom = L.Vmom;           % dipole orientations

% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
    E.A{i} = log(A{i} + eps);                      % forward
    V.A{i} = A{i}/2;                               % backward
    Q      = Q | A{i};                             % and lateral connections
end

for i = 1:length(B)
    E.B{i} = 0*B{i};                               % input-dependent scaling
    V.B{i} = B{i}/2;
    Q      = Q | B{i};
end
E.C        = log(C + eps);                         % where inputs enter
V.C        = C/2;

% set delay (enforcing symmetric delays)
%--------------------------------------------------------------------------
E.D        = sparse(n,n);
V.D        = Q/16;

% set stimulus parameters
%--------------------------------------------------------------------------
E.R        = [0 0];        V.R = [1 1]/16;         % input [Gamma] parameters
E.N        = sparse(m,1);  V.N = ones(m,1);        % DCT coefficients
E.U        = U;            V.U = 0;                % trial duration (s)

% vectorize
%--------------------------------------------------------------------------
pE         = E;
pV         = spm_vec(V);
pC         = diag(sparse(pV));
warning on



% prior momments if two arguments
%--------------------------------------------------------------------------
if nargout == 2, varargout{1} = pE; varargout{2} = pC; return, end


% Model specification
%==========================================================================
clear global M
global M

M.IS     = 'spm_int_U';
M.f      = 'spm_fx_erp';
M.g      = 'spm_gx_erp';
M.x      = sparse(n*9 + 1,1);
M.pE     = pE;
M.pC     = pC;
M.m      = length(B);
M.n      = length(M.x);
M.l      = n;
M.Nareas = n;
M.Spatial_type = 3;

if nargout == 1, varargout{1} = M; return, end

% compute impulse response
%--------------------------------------------------------------------------
clear U
M.pE.K  = 1;
M.onset = 128;
N       = 128;
U.u     = sparse(N,M.m);
U.dt    = 8/1000;
y       = feval(M.IS,M.pE,M,U);
plot([1:N]*U.dt*1000,y)


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



