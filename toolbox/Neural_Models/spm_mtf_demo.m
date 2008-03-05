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
k     = k(f(k) > 2 & f(k) < 64);
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
B     = {};
C     = sparse(n,1,1,n,1);
 
% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(A,B,C);
 
% create LFP model
%--------------------------------------------------------------------------
M.IS  = 'spm_lfp_mtf';
M.FS  = 'spm_lfp_sqrt';
M.f   = 'spm_fx_lfp';
M.g   = 'spm_gx_erp';
M.x   = sparse(n,13);
M.n   = length(M.x(:));
M.pE  = pE;
M.pC  = pC;
M.m   = length(B);
M.l   = 1;
M.Hz  = f;

 
% inversion (in frequency space)
%==========================================================================
 
% data
%--------------------------------------------------------------------------
y      = spm_cond_units(y)*32;
Y.y    = {y}; 
Y.Q    = {spm_Q(1/2,length(f),1)};
 
% invert
%--------------------------------------------------------------------------
[Ep,Cp,S,F] = spm_nlsi_GN(M,[],Y);
 
% plot spectral density (after removing mean in log-space)
%--------------------------------------------------------------------------
[G w]  = spm_lfp_mtf(Ep,M);
 
subplot(2,1,1)
plot(w,G{1},w,y,':')
xlabel('frequency (Hz)')
xlabel('Power')
legend({'predicted','observed'})
axis square
grid on
