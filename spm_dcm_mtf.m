function [S,K,s,w,t,dfdx] = spm_dcm_mtf(P,M,U)
% computes transfer functions using the system's eigenspectrum
% FORMAT [S,K,s,w,t,dfdx] = spm_dcm_mtf(P,M,[U])
%
% P - model parameters
% M - model (with flow M.f and expansion point M.x and M.u)
% U - induces expansion around steady state (from spm_dcm_neural_x(P,M))
%
% S    - modulation transfer functions (complex)
% K    - Volterra kernels (real)
% s    - eigenspectrum (complex)
% w    - frequencies (Hz) = M.Hz
% t    - time (seconds)   = M.pst
% dfdx - Jacobian
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
% $Id: spm_dcm_mtf.m 6856 2016-08-10 17:55:05Z karl $


% get local linear approximation
%==========================================================================

% solve for steady-state - if required
%--------------------------------------------------------------------------
if nargin > 2
    M.x   = spm_dcm_neural_x(P,M);
end

% check expansion points
%--------------------------------------------------------------------------
try, M.x; catch, M.x = spm_dcm_x_neural(P,M.dipfit.model); end
try, M.u; catch, M.u = sparse(M.m,1);                      end

% frequencies and peristimulus time of interest
%--------------------------------------------------------------------------
w      = (1:64)';
t      = (0:64 - 1)'/64;
try, w = M.Hz(:);       end
try, t = M.pst(:);      end
try, t = M.dt*(1:M.N)'; end


% delay operator - if not specified already
%--------------------------------------------------------------------------
if isfield(M,'D')
    
    dfdx = spm_diff(M.f,M.x,M.u,P,M,1);
    dfdu = spm_diff(M.f,M.x,M.u,P,M,2);
    D    = M.D;
    
else
    
    if nargout(M.f) == 4
        [f,dfdx,D,dfdu] = feval(M.f,M.x,M.u,P,M);
        
    elseif nargout(M.f) == 3
        [f,dfdx,D]      = feval(M.f,M.x,M.u,P,M);
        dfdu            = spm_diff(M.f,M.x,M.u,P,M,2);
        
    elseif nargout(M.f) == 2
        [f,dfdx]        = feval(M.f,M.x,M.u,P,M);
        dfdu            = spm_diff(M.f,M.x,M.u,P,M,2);
        D               = 1;
    else
        dfdx            = spm_diff(M.f,M.x,M.u,P,M,1);
        dfdu            = spm_diff(M.f,M.x,M.u,P,M,2);
        D               = 1;
    end
end

% Jacobian and eigenspectrum
%==========================================================================
if nargout(M.g) == 2
    [g,dgdx] = feval(M.g,M.x,M.u,P,M);
else
    dgdx     = spm_diff(M.g,M.x,M.u,P,M,1);
end
dfdx  = D*dfdx;
dfdu  = D*dfdu;
try
   [v,s] = eig(full(dfdx),'nobalance');
catch
   v  = eye(size(dfdx));
   s  = NaN(size(dfdx));
end
s     = diag(s);


% condition unstable eigenmodes
%--------------------------------------------------------------------------
if max(w) > 1
    s = 1j*imag(s) + real(s) - exp(real(s));
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
            Sk       = 1./(1j*2*pi*w - s(k));
            S(:,i,j) = S(:,i,j) + dgdv(i,k)*dvdu(k,j)*Sk;
            
            % kernels
            %--------------------------------------------------------------          
            if nargout > 1
                Kk       = exp(s(k)*t);
                K(:,i,j) = K(:,i,j) + real(dgdv(i,k)*dvdu(k,j)*Kk);
            end
            
        end
    end
end



return

% NOTES: internal consistency with explicit Fourier transform of kernels
%==========================================================================

% augment and bi-linearise (with intrinsic delays)
%--------------------------------------------------------------------------
M.D       = D;
[M0,M1,L] = spm_bireduce(M,P);

% project onto spatial modes
%--------------------------------------------------------------------------
try,    L = M.U'*L; end

% kernels
%--------------------------------------------------------------------------
N         = length(t);
dt        = (t(2) - t(1));
[K0,K1]   = spm_kernels(M0,M1,L,N,dt);

% Transfer functions (FFT of kernel)
%--------------------------------------------------------------------------
S1        = fft(K1)*dt;
w1        = ((1:N) - 1)/(N*dt);
j         = w1 < max(w);
S1        = S1(j,1,1);
w1        = w1(j);

subplot(2,2,1), plot(t,K(:,1,1),t,K1(:,1,1));
title('kernels','fontsize',16)
xlabel('peristimulus time')


subplot(2,2,2), plot(w,real(S(:,1,1)),    w1,real(S1(:,1,1))); hold on
subplot(2,2,2), plot(w,imag(S(:,1,1)),':',w1,imag(S1(:,1,1)),':'); hold off
title('transfer functions','fontsize',16)
xlabel('frequency');








