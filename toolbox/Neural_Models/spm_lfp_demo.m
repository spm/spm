% Demo routine for local field potential models
%==========================================================================

clear global
clear

% number of regions in coupled map lattice
%--------------------------------------------------------------------------
n     = 1;

% specifc network (connections)
%--------------------------------------------------------------------------
if n > 1
A{1}  = diag(ones(n - 1,1),-1);
else
    A{1} = 0;
end
A{2}  = A{1}';
A{3}  = sparse(n,n);
B{1}  = sparse(n,n);
C     = sparse(1,1,256,n,1);

% mixture of regions subtending LFP(EEG)
%--------------------------------------------------------------------------
L     = ones(1,n);

% mixture of states subtending LFP(EEG)
%--------------------------------------------------------------------------
H      = sparse(13,2);
H(9,1) = 1;               % pyramidal depolorization
H(1,2) = 1;               % stellate  depolorization

% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(A,B,C,L,H);

% create LFP model
%--------------------------------------------------------------------------
[l n] = size(L);

M.f   = 'spm_fx_lfp';
M.g   = 'spm_gx_lfp';
M.x   = sparse(n,13);
M.pE  = pE;
M.pC  = pC;
M.m   = size(C,2);
M.n   = length(M.x(:));
M.l   = size(H,2);
M.IS  = 'spm_int_J';

% create BOLD model
%--------------------------------------------------------------------------
[hE,hC] = spm_hdm_priors(1,5);
hE(end) = 1;

% model
%--------------------------------------------------------------------------
clear H
H.f     = 'spm_fx_hdm';
H.g     = 'spm_gx_hdm';
H.x     = [0 0 0 0]';
H.pE    = hE;    
H.pC    = hC;
H.m     = 1;
H.n     = 4;
H.l     = 1;

% for Andre and Jean
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


% Integrate system to see response
%--------------------------------------------------------------------------
N     = 2048;
U.dt  = 8/1000;
U.u   = 32*(sparse(128:512,1,1,N,M.m) + randn(N,M.m)/16);
t     = [1:N]*U.dt;
LFP   = spm_int_J(pE,M,U);

% input
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(t,U.u)
title('Exogeneous input')
axis square
xlabel('time (s)')

% LFP
%--------------------------------------------------------------------------
subplot(2,2,2)
title('LFP response')
plot(t,LFP)
axis square
xlabel('time (s)')

% time-frequency
%--------------------------------------------------------------------------
W     = 512;
w     = 4:1/4:32;
cpW   = w*W*U.dt;
subplot(2,2,3)
imagesc(t,w,abs(spm_wft(LFP(:,1),cpW,W)).^2);
title('time-frequency response')
axis square xy
xlabel('time (s)')
ylabel('Hz')

% Use response to drive a hemodynamic model
%--------------------------------------------------------------------------
sigmoid = inline('1./(1 + exp(-2*(x - 1))) - 1./(1 + exp(2))','x');
U.u   = sigmoid(U.u);
BOLD  = spm_int_J(hE,H,U);

subplot(2,2,4)
plot(t,BOLD)
title('BOLD response')
axis square
xlabel('time (s)')


% Stability analysis (over excitatory and inhibitory time constants)
%--------------------------------------------------------------------------
p     = [-16:16]/8;
np    = length(p);
HZ    = sparse(np,np);
for i = 1:np
    for j = 1:np
        P        = M.pE;
        P.T(:,1) = P.T(:,1) + p(i);
        P.T(:,2) = P.T(:,2) + p(j);
        S        = eig(full(spm_bireduce(M,P)));
        LE(i,j)  = max(real(S));
        S        = S(abs(imag(S)) > 2*2*pi & abs(imag(S)) < 64*2*pi);
        try
            [k l]    = max(real(S));
            HZ(i,j)  = imag(S(l))/(2*pi);
        end
    end
end

p1  = 4*exp(p);
p2  = 16*exp(p);

subplot(2,2,1)
surf(p1,p2,LE')
shading interp
axis square
ylabel('inbitory time constant (ms)')
xlabel('excitatory time constant (ms)')
title('stability')

subplot(2,2,2)
contour(p2,p1,LE,[0 0])
axis square
xlabel('inbitory time constant')
ylabel('excitatory time constant')
title('stability')

subplot(2,2,3)
surf(p1,p2,HZ')
shading interp
axis square
ylabel('inbitory time constant')
xlabel('excitatory time constant')
title('Frequency')

subplot(2,2,4)
contour(p2,p1,full(HZ),[8 12 16 30],'k');
axis square
xlabel('inbitory time constant')
ylabel('excitatory time constant')
title('Frequency')


% transfer functions
%==========================================================================

% compute transfer functions (switch off noise)
%--------------------------------------------------------------------------
pE    = M.pE;
pE.a  = -32;
pE.b  = -32;
[G w] = spm_lfp_mtf(pE,M);

subplot(2,1,1)
plot(w,G)
axis square
xlabel('frequency {Hz}')


% compute transfer functions for different inhibitory time constants
%--------------------------------------------------------------------------
p     = log([1:64]/32);
for i = 1:length(p)
    Pi      = pE;
    Pi.T(2) = Pi.T(2) + p(i);
    GW(:,i) = spm_lfp_mtf(Pi,M);
end

subplot(2,1,1)
imagesc(16*exp(p),w,GW)
ylabel('Frequency')
xlabel('Inhibitory time constant (ms)')

subplot(2,1,2)
plot(w,GW)
xlabel('Frequency')
ylabel('g(w)')



% Integrate system to see Transient response
%==========================================================================
pE.T(2) = log(2);
N     = 1024;
U.dt  = 1/1000;
U.u   = 8*(exp(-([1:N]' - N/4).^2/(2*32^2)) + rand(N,1)/4);
t     = [1:N]*U.dt;
LFP   = spm_int_J(pE,M,U);

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
W     = 512;
w     = 4:1/4:32;
cpW   = w*W*U.dt;
subplot(2,2,3)
imagesc(t*1000,w,abs(spm_wft(LFP(:,1),cpW,W)).^2);
title('time-frequency response')
axis square xy
xlabel('time (ms)')





