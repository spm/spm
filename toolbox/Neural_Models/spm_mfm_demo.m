% Demo routine for mean-field models
%==========================================================================

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
clear M
n     = 1;

% extrinsic network connections
%--------------------------------------------------------------------------
if n > 1
A{1}  = diag(ones(n - 1,1),-1);
else
    A{1} = 0;
end
A{2}  = A{1}';
A{3}  = sparse(n,n);
B{1}  = sparse(n,n);
C     = sparse(1,1,1,n,1);


% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_mfm_priors(A,B,C);

% create LFP model
%--------------------------------------------------------------------------
M.f   = 'spm_fx_mfm';
M.x   = spm_x_mfm(pE);
M.pE  = pE;
M.pC  = pC;
M.m   = size(C,2);
M.n   = length(spm_vec(M.x));
M.l   = size(L,1);

% solve for fixed point 
%--------------------------------------------------------------------------
U.u   = sparse(128,1);
U.dt  = 4/1000;
x     = spm_int_ode(pE,M,U);
x     = spm_unvec(x(end,:),M.x);
M.x   = x;


% Integrate system to see Transient response
%==========================================================================
%M.g   = 'spm_gx_erp';
N     = 256;
U.u   = 8*(exp(-([1:N]' - N/4).^2/(2*4^2)) + rand(N,1)/exp(32));
t     = [1:N]*U.dt;
LFP   = spm_int_B(pE,M,U);

% LFP
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(t*1000,U.u)
title('input')
axis square
xlabel('time (ms)')

% LFP
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(t*1000,LFP)
title('depolarization')
axis square
xlabel('time (ms)')

% time-frequency
%--------------------------------------------------------------------------
W     = 256;
w     = 4:1/4:32;
cpW   = w*W*U.dt;
subplot(2,2,3)
imagesc(t*1000,w,abs(spm_wft(LFP(:,7),cpW,W)).^2);
title('time-frequency response')
axis square xy
xlabel('time (ms)')


return


% Kernels
%==========================================================================

% augment and bi-linearise
%--------------------------------------------------------------------------
[M0,M1,L1,L2] = spm_bireduce(M,M.pE);

% compute kernels (over 64 ms)
%--------------------------------------------------------------------------
N          = 64;
dt         = 1/1000;
t          = [1:N]*dt*1000;
[K0,K1,K2] = spm_kernels(M0,M1,L1,L2,N,dt);

subplot(2,1,1)
plot(t,K1)
title('1st-order Volterra kernel')
axis square
xlabel('time (ms)')

subplot(2,1,2)
imagesc(t,t,K2(1:64,1:64,1))
title('2nd-order Volterra kernel')
axis square
xlabel('time (ms)')






