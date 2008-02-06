% Demo for models of induced frequency responses and nonlinear coupling
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
    A{1}  = 0;
end
A{2}  = A{1}';
A{3}  = sparse(n,n);
B{1}  = sparse(n,n);
C     = sparse(1,1,1,n,1);
 
% mixture of regions subtending LFP(EEG)
%--------------------------------------------------------------------------
L     = ones(1,n);
 
% mixture of states subtending LFP(EEG)
%--------------------------------------------------------------------------
H     = sparse(9,1,1,13,1);
 
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
M.m   = length(B);
M.n   = size(M.x,1);
M.l   = l;
 
% Integrate system to see response
%--------------------------------------------------------------------------
N     = 1024;                         % number of samples
U.dt  = 8/1000;                       % sampling interval
pst   = [1:N]*U.dt;                   % peristimulus time
t     = 32;                           % sample window for WFT
cpt   = 1:1/8:16;                     % cycles per window         
w     = cpt./(t*U.dt);                % Hz
 
% input 
%==========================================================================
U.u   = sparse(128:512,1,1,N,1)*64 + randn(N,1)*4;        % noisy burst
U.u   = sparse(128:512,1,1,N,1).*sin(2*pi*16*pst(:));     % pure Hz - low
U.u   = sparse(128:512,1,1,N,1).*sin(2*pi*16*pst(:))*128; % pure Hz - high
 
% response
%--------------------------------------------------------------------------
LFP   = spm_int_ode(pE,M,U);
 
% display
%==========================================================================
 
% input - time
%--------------------------------------------------------------------------
subplot(2,2,1)
plot(pst,U.u)
axis square tight
xlabel('time (s)')
ylabel('activity')
 
% Input - time-frequency
%--------------------------------------------------------------------------
subplot(2,2,3)
imagesc(pst,w,sqrt(abs(spm_wft(U.u,cpt,t))));
axis square xy
xlabel('time (s)')
ylabel('frequency')
 
% LFP - time
%--------------------------------------------------------------------------
subplot(2,2,2)
plot(pst,LFP)
axis square tight
xlabel('time (s)')
ylabel('activity')
 
% LFP - time-frequency
%--------------------------------------------------------------------------
subplot(2,2,4)
imagesc(pst,w,sqrt(abs(spm_wft(LFP,cpt,t))));
axis square xy
xlabel('time (s)')
ylabel('frequency')
