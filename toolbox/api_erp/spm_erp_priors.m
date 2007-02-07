function [varargout] = spm_erp_priors(A,B,C,dipfit,dur,ons)
% prior moments for a neural mass model of erps
% FORMAT [pE,pC] = spm_erp_priors(A,B,C,dipfit,dur,ons)
% FORMAT [M]     = spm_erp_priors(A,B,C,dipfit,dur,ons)
% FORMAT           spm_erp_priors(A,B,C,dipfit,dur,ons)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
% dipfit       - prior forward model structure
% dur          - trial duration (sec) [default = 1 s]
% ons          - onset (nsec) [default = 60 ms]
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
% stimulus and noise parameters
%--------------------------------------------------------------------------
%    pE.R - onset and dispersion (Gamma paramters)
%    pE.M - magnitude
%    pE.N - background fluctuations
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

% default: a single source model
%--------------------------------------------------------------------------
if nargin < 3
    A   = {0 0 0};
    B   = {};
    C   = 1;
end
if nargin < 4
    dipfit.L.pos = sparse(3,0);
end
if nargin < 5
    dur = 1;
end
if nargin < 6
    ons = 60;
end
try
    dipfit.type;
catch
    dipfit.type = 'ECD';
end

% disable log zero warning
%--------------------------------------------------------------------------
warning off
n     = size(C,1);                                 % number of sources
n1    = ones(n,1);

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(n1);        V.T = n1/8;                % time constants
E.H   = log(n1);        V.H = n1/8;                % synaptic density

% set intrinic [excitatory] time constants
%--------------------------------------------------------------------------
E.S   = [0 0];          V.S = [1 1]/8;             % dispersion & threshold

% paramters for electromagnetic forward model
%--------------------------------------------------------------------------
if strcmp(dipfit.type,'ECD')
    
    E.Lpos = dipfit.L.pos;  V.Lpos =  16*ones(3,n);    % dipole positions
    E.Lmom = sparse(3,n);   V.Lmom = 256*ones(3,n);    % dipole orientations

else
    m      = dipfit.Nm;
    E.L    = sparse(m,n);   V.Lpos =  16*ones(m,n);    % dipole modes
end

% set extrinsic connectivity
%--------------------------------------------------------------------------
Q     = sparse(n,n);
for i = 1:length(A)
    E.A{i} = log(A{i} + eps);                      % forward
    V.A{i} = A{i};                                 % backward
    Q      = Q | A{i};                             % and lateral connections
end

for i = 1:length(B)
    E.B{i} = 0*B{i};                               % input-dependent scaling
    V.B{i} = B{i};
    Q      = Q | B{i};
end
E.C        = log(C + eps);                         % where inputs enter
V.C        = C;

% set delay (enforcing symmetric delays)
%--------------------------------------------------------------------------
E.D        = sparse(n,n);
V.D        = Q/8;

% set initial states (Spiny stellate depolarisarion))
%--------------------------------------------------------------------------
% E.x0     = sparse(n,1);
% V.x0     = n1;

% set stimulus parameters onset, dispersion and duration
%--------------------------------------------------------------------------
E.dur      = dur;                V.dur = 0;            % trial duration  (ms)
E.ons      = ons;                V.ons = 0*ons;        % stimulus onsets (ms)
u          = length(ons);
E.R        = sparse(u,2);        V.R   = ones(u,2)/8;  % input [Gamma] parameters
E.M        = sparse(1,1,1,u,1);  V.M   = ones(u,1);    % magnitude



% background fluctuations
%--------------------------------------------------------------------------
E.N        = ones(n,1)*[0 0 12]; V.N   = ones(n,3)*8;  % amplitude and Hz

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
M.ns     = 128;
M.dipfit.type = 'LFP';

if nargout == 1, varargout{1} = M; return, end

% compute impulse response
%--------------------------------------------------------------------------
clear U
U.dt    = 4/1000;
y       = feval(M.IS,M.pE,M,U);
plot([1:M.ns]*U.dt*1000,y)


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
axis square



