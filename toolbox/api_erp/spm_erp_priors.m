function [pE,pC] = spm_erp_priors(A,B,C,L,U)
% prior moments for a neural mass model of erps
% FORMAT [pE,pC] = spm_erp_priors(A,B,C,L,U)
%
% A{3},B{m},C    - binary constraints on extrinsic connectivity
% L              - lead field matrix
% U              - trial duration (sec)      [default = 1 s] 
%
% pE - prior expectation
% pC - prior covariances
%
% Because priors are specified under log normal assumptions, most
% parameters are simply scaling coefficients with a prior expectation
% and variance of one.  After log transform this renders pE = 0 and
% pC = 1;  The prior expectations of what they scale are specified in
% spm_erp_fx
%___________________________________________________________________________
%
% David O, Friston KJ (2003) A neural mass model for MEG/EEG: coupling and
% neuronal dynamics. NeuroImage 20: 1743-1755
%___________________________________________________________________________
% %W% Karl Friston %E%

% impulse response for a single area if no input arguments
%---------------------------------------------------------------------------
if nargin == 0, spm_erp,        return, end
if nargin == 3, spm_erp(A,B,C), return, end

% default duration
%---------------------------------------------------------------------------
if nargin <  5, U = 1; end
warning off

% disable log zero warning
%---------------------------------------------------------------------------
n     = size(C,1);                                 % number of sources
m     = 8;                                         % and noise components
n1    = ones(n,1);

% set intrinic [excitatory] time constants
%---------------------------------------------------------------------------
E.T   = log(n1);        V.T = n1/16;               % time constants
E.H   = log(n1);        V.H = n1/16;               % synaptic density

% set observer parameters
%---------------------------------------------------------------------------
if ~isstruct(L)
    % static leadfield
    E.K   = n1*0;           V.K = sparse(n1);          % pyramidal CSD
    E.L   = L;              V.L = sparse(0*L);         % lead field
else
   % parameterised leadfield based on equivalent current dipoles 
   % prior distributions given by user
   E.K = L.K;                   V.K = L.VK;
   E.Lpos = L.pos;              V.Lpos = L.Vpos;         % dipole positions
   E.Lmom = L.mom;              V.Lmom = L.Vmom;         % dipole orientations
end

% set extrinsic connectivity
%---------------------------------------------------------------------------
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
%---------------------------------------------------------------------------
E.D        = sparse(n,n);
V.D        = Q/16;

% set stimulus parameters
%---------------------------------------------------------------------------
E.R        = [0 0];        V.R = [1 1]/16;         % input [Gamma] parameters
E.N        = sparse(m,1);  V.N = ones(m,1);        % DCT coefficients
E.U        = U;            V.U = 0;                % trial duration (s)

% vectorize
%---------------------------------------------------------------------------
pE         = spm_vec(E);
pV         = spm_vec(V);
pC         = diag(sparse(pV));
warning on

return


% Impulse response
%===========================================================================
function spm_erp(A,B,C)

% model specification
%---------------------------------------------------------------------------
if ~nargin                                                 % a single source
    A   = {0 0 0};
    B   = {};
    C   = 1;
end
n       = size(C,1);
[qE,qC] = spm_erp_priors(A,B,C,[]);
IS    = 'spm_int_U';
M.f   = 'spm_fx_erp';
M.g   = 'spm_gx_erp';
M.x   = sparse(n*9 + 1,1);
M.pE  = qE;
M.pC  = qC;
M.m   = length(B);
M.n   = length(M.x);
M.l   = n*2;
M.IS  = IS;

% compute impulse response
%---------------------------------------------------------------------------
N     = 64;
U.u   = sparse(N,M.m);
U.dt  = 1/N/2;
y     = feval(IS,qE,M,U);
plot([1:N]*U.dt*1000,y)


return

% demo for log-normal pdf
%---------------------------------------------------------------------------
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



