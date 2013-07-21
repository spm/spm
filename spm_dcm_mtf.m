function [S,K,s,w,t] = spm_dcm_mtf(P,M,U)
% computes transfer functions using the system's eigenspectrum
% FORMAT [S,K,s,w,t] = spm_dcm_mtf(P,M,[U])
%
% P - model parameters
% M - model (with flow M.f and expansion point M.x and M.u)
% U - induces expansion around steady state (from spm_dcm_neural_x(P,M))
%
% S - directed transfer functions (complex)
% K - directed kernels (real)
% s - eigenspectrum (complex)
% w - frequencies (Hz) = M.Hz
% t - time (seconds)   = M.pst
%
% This routine uses the eigensolution of a dynamical systems Jacobian to
% complete the first-order Volterra terminals and transfer functions  in
% peristimulus and frequency space respectively.  The advantage of using
% the-solution is that unstable modes (eigenvectors of the Jacobian) can be
% conditioned (suppressed). Furthermore, this provides for a
% computationally efficient and transparent evaluation of the transfer
% functions that draws on linear signal processing theory in frequency
% space.
%__________________________________________________________________________
% Copyright (C) 2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_mtf.m 5588 2013-07-21 20:59:39Z karl $


% get local linear approximation
%==========================================================================

% solve for steady-state - if required
%--------------------------------------------------------------------------
if nargin > 2
    M.x   = spm_dcm_neural_x(P,M);
end

% check expansion points
%--------------------------------------------------------------------------
try, M.x; catch, M.x = sparse(M.n,1); end
try, M.u; catch, M.u = sparse(M.l,1); end

% frequencies and peristimulus time of interest
%--------------------------------------------------------------------------
w      = (1:64)';
t      = (0:64 - 1)'/64;
try, w = M.Hz(:);       end
try, t = M.pst(:);      end
try, t = M.dt*(1:M.N)'; end


% delay operator - if not specified already
%--------------------------------------------------------------------------
if nargout(M.f) == 3
    [f,dfdx,D] = feval(M.f,M.x,M.u,P,M);
    
elseif nargout(M.f) == 2
    [f,dfdx]   = feval(M.f,M.x,M.u,P,M);
    D          = 1;
else
    dfdx       = spm_diff(M.f,M.x,M.u,P,M,1);
    D          = 1;
end

% Jacobian and eigenspectrum
%==========================================================================
dfdu  = spm_diff(M.f,M.x,M.u,P,M,2);
dgdx  = spm_diff(M.g,M.x,M.u,P,M,1);
dfdx  = D*dfdx;
dfdu  = D*dfdu;

[v,s] = eig(full(dfdx),'nobalance');
s     = diag(s);


% condition unstable eigenmodes
%--------------------------------------------------------------------------
if max(w) > 1
    s = 1j*imag(s) + min(real(s),-4);
else
    s = 1j*imag(s) + min(real(s),-1/32);
end


% Transfer functions
%==========================================================================

% transfer functions (FFT of kernel)
%--------------------------------------------------------------------------
nw    = size(w,1);            % number of frequencies
nt    = size(t,1);            % number of time bins
ng    = size(dgdx,1);         % number of outputs
nu    = size(dfdu,2);         % number of inputs
nk    = size(v,2);            % number of modes
S     = zeros(nw,ng,nu);
K     = zeros(nt,ng,nu);

% derivatives over modes
%--------------------------------------------------------------------------
dgdv  = dgdx*v;
dvdu  = pinv(v)*dfdu;

for j = 1:nu
    for i = 1:ng
        for k = 1:nk
            
            % transfer functions (FFT of kernel)
            %--------------------------------------------------------------
            Sk       = 1./(1j*(2*pi*w - imag(s(k))) - real(s(k)));
            S(:,i,j) = S(:,i,j) + dgdv(i,k)*dvdu(k,j)*Sk;
            
            if nargout > 1
                
                % kernels
                %----------------------------------------------------------
                Kk       = exp(s(k)*t);
                K(:,i,j) = K(:,i,j) + real(dgdv(i,k)*dvdu(k,j)*Kk);
                
            end
        end
    end
end

