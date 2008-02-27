function [varargout] = spm_erp_priors(A,B,C,dipfit,u)
% prior moments for a neural-mass model of ERPs
% FORMAT [pE,gE,pC,gC] = spm_erp_priors(A,B,C,dipfit)
% FORMAT [pE,pC] = spm_erp_priors(A,B,C,dipfit)
% FORMAT [M]     = spm_erp_priors(A,B,C,dipfit)
% FORMAT           spm_erp_priors(A,B,C,dipfit)
%
% A{3},B{m},C  - binary constraints on extrinsic connections
% dipfit       - prior forward model structure
% u            - number of neuronal inputs
%
% pE - prior expectation - f(x,u,P,M)
% gE - prior expectation - g(x,u,G,M)
%
% synaptic parameters
%--------------------------------------------------------------------------
%    pE.T - syaptic time constants
%    pE.H - syaptic densities
%    pE.S - activation function parameters
%
% spatial parameters
%--------------------------------------------------------------------------
%    gE.Lpos - position                   - ECD
%    gE.L    - orientation
%
% or gE.L    - coefficients of local modes - Imaging
%    gE.L    - gain of electrodes         - LFP
%
% source contribution to ECD
%--------------------------------------------------------------------------
%    gE.J    - contribution
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
%    pE.R - onset and dispersion
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
% $Id: spm_erp_priors.m 1174 2008-02-27 20:22:30Z karl $
 
% default: a single source model
%--------------------------------------------------------------------------
if nargin < 3
    A   = {0 0 0};
    B   = {};
    C   = 1;
end
if nargin < 5, u = 1; end
 
% disable log zero warning
%--------------------------------------------------------------------------
warning off
n     = size(C,1);                                   % number of sources
n1    = ones(n,1);
 
% paramters for electromagnetic forward model
%==========================================================================
try,   type = dipfit.type; catch, type = 'LFP'; end
switch type
    
    case{'ECD (EEG)','ECD (MEG)'}
    %----------------------------------------------------------------------
    G.Lpos = dipfit.L.pos;  U.Lpos =      0*ones(3,n);    % positions
    G.L    = sparse(3,n);   U.L    = exp(8)*ones(3,n);    % orientations
 
    case{'Imaging'}
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
 
% set intrinsic [excitatory] time constants
%--------------------------------------------------------------------------
E.T   = log(n1);        V.T = n1/8;                % time constants
E.H   = log(n1);        V.H = n1/32;               % synaptic density
 
% set parameter of activation function
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
V.D        = Q/16;
 
% set stimulus parameters: onset and dispersion
%--------------------------------------------------------------------------
E.R        = sparse(u,2);  V.R   = ones(u,1)*[1/16 1/16];
warning on
 
% prior moments if two arguments
%==========================================================================
if nargout == 4
    varargout{1} = E;
    varargout{2} = G;
    varargout{3} = diag(sparse(spm_vec(V)));
    varargout{4} = diag(sparse(spm_vec(U)));
    return
else
    switch ype
 
        case{'ECD (EEG)','ECD (MEG)'}
            %--------------------------------------------------------------
            E.Lpos  = G.Lpos;
            E.L     = G.L   ;
            C       = diag(sparse(spm_vec(V,U)));
 
        case{'Imaging','LFP'}
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
 
 
% Model specification
%==========================================================================
M.IS          = 'spm_int_U';
M.f           = 'spm_fx_erp';
M.g           = 'spm_gx_erp';
M.x           = sparse(n*9 + 1,1);
M.pE          = E;
M.pC          = C;
M.m           = length(B);
M.n           = length(M.x);
M.l           = n;
M.ns          = 128;
M.ons         = 80;
M.dipfit.type = 'LFP';
 
if nargout == 1, varargout{1} = M; return, end
 
% compute impulse response
%--------------------------------------------------------------------------
clear U
U.dt    = 1/1000;
M.ns    = 256;
y       = feval(M.IS,M.pE,M,U);
plot([1:M.ns]*U.dt*1000,y)
 
 
return
 
% demo for log-normal pdf
%==========================================================================
x    = [1:64]/16;
for i = [2 16 128]
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
