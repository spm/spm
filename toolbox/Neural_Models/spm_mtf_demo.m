% Demo routine for inverting local field potential models
%==========================================================================
 
clear global
clear
 
% emprical data - sort and decimate
%--------------------------------------------------------------------------
load  'ten_minute_average_control.mat';
y     = G_control(:);
f     = f_Control(:);
[f i] = sort(f);
y     = y(i);
 
for i = 1:64
    [d j] = min(abs(f - i));
    k(i)  = j;
end
f     = f(k);                  % frequency
y     = y(k);                  % power
 
 
% specify model
%==========================================================================
 
% number of regions in coupled map lattice
%--------------------------------------------------------------------------
n     = 1;
 
% specifc network (connections)
%--------------------------------------------------------------------------
A{1}  = triu(ones(n,n),1);
A{2}  = sparse(n,n);
A{3}  = sparse(n,n);
B{1}  = sparse(n,n);
C     = sparse(n,1,1,n,1);
 
% mixture of regions subtending LFP(EEG)
%--------------------------------------------------------------------------
L     = sparse(1,1,1,1,n);
 
% mixture of states subtending LFP(EEG)
%--------------------------------------------------------------------------
H     = sparse(9,1,1,13,1);
 
% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(A,B,C,L,H);
 
% create LFP model
%--------------------------------------------------------------------------
[l n] = size(L);
 
M.IS  = 'spm_lfp_mtf';
M.FS  = 'log';
M.f   = 'spm_fx_lfp';
M.g   = 'spm_gx_lfp';
M.x   = sparse(n,13);
M.n   = length(M.x(:));
M.pE  = pE;
M.pC  = pC;
M.m   = length(B);
M.l   = l;
M.Hz  = f;
 
% Gain
%--------------------------------------------------------------------------
[G w]  = spm_lfp_mtf(pE,M);
 
% inversion (in frequency space)
%==========================================================================
 
% data
%--------------------------------------------------------------------------
Y.y    = y*norm(G)/norm(y); 
Y.X0   = ones(length(y),1);
Y.Q    = {spm_Q(2/3,length(f),1)};
 
% invert
%--------------------------------------------------------------------------
[Ep,Cp,S,F] = spm_nlsi_GN(M,[],Y);
 
% plot spectral density (after removing mean in log-space)
%--------------------------------------------------------------------------
[G w]  = spm_lfp_mtf(Ep,M);
G      = exp(log(G)   - mean(log(G)));
y      = exp(log(Y.y) - mean(log(Y.y)));
 
subplot(2,1,1)
plot(w,G,w,y,':')
xlabel('frequency (Hz)')
xlabel('Power')
legend({'predicted','observed'})
axis square
grid on
 
 
% Plot posterior Estimates
%--------------------------------------------------------------------------
PRIOR.R = exp(pE.R).*[2 1];
PRIOR.T = exp(pE.T).*[4 16];
PRIOR.H = exp(pE.H).*[4 32];
PRIOR.G = exp(pE.G).*[128 128 64 64 16];
PRIOR.I = exp(pE.I).*[2];
 
MAP.R   = exp(Ep.R).*[2 1];
MAP.T   = exp(Ep.T).*[4 16];
MAP.H   = exp(Ep.H).*[4 32];
MAP.G   = exp(Ep.G).*[128 128 64 64 16];
MAP.I   = exp(Ep.I).*[2];
 
subplot(2,1,2)
bar([spm_vec(MAP) spm_vec(PRIOR)]);
title('MAP Estimates');
xlabel('Parameters');
ylabel('Conditional mean');
legend({'posterior','prior'})
set(gca,'Xticklabel',{'gain','shift','T-ex','T-in','H-ex','H-in','g1','g2','g3','g4','g5','delay'})
