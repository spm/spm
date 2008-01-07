% Demo routine for inverting local field potential models using
% cross-spectral density summaries of steady-state dynamics
%==========================================================================
 
clear global
clear
 
% specify model
%==========================================================================
 
% number of sources and LFP channels (usually the same)
%--------------------------------------------------------------------------
n     = 2; % number of sources
nc    = 2; % number of channels
 
% specify network (connections)
%--------------------------------------------------------------------------
A{1}  = tril(ones(n,n),-1);                % a forward connection
A{2}  = triu(ones(n,n),+1);                % a backward connection
A{3}  = sparse(n,n);                       % lateral connections
for i = 1:n
    B{i} = sparse(n,n);                    % input-specific modulation
end                                        % not used in this context
C     = speye(n,n);                        % sources receiving innovations
 
 
% mixture of states subtending LFP(EEG) from each source
%--------------------------------------------------------------------------
H     = sparse(9,1,1,13,1);
 
% mixture of sources subtending LFP(EEG); c.f., a lead field
%--------------------------------------------------------------------------
L     = speye(nc,n);
 
% get priors
%--------------------------------------------------------------------------
[pE,pC] = spm_lfp_priors(A,B,C,L,H);
 
% create LFP model
%--------------------------------------------------------------------------
M.IS  = 'spm_lfp_mtf';
M.FS  = 'spm_lfp_log';
M.f   = 'spm_fx_lfp';
M.g   = 'spm_gx_lfp';
M.x   = sparse(n,13);
M.n   = n*13;
M.pE  = pE;
M.pC  = pC;
M.m   = n;
M.l   = nc;
M.Hz  = [1:64]';                          % frequencies
 
 
% simulate spectral data directly
%==========================================================================
P           = pE;
P.A{1}(2,1) = 1;                          % strong forward connections
CSD         = spm_lfp_mtf(P,M);
 
% or generate data
%==========================================================================
 
% Integrate with pink noise process
%--------------------------------------------------------------------------
N     = 512;
U.dt  = 8/1000;
U.u   = randn(N,M.m)/16;
U.u   = sqrt(spm_Q(1/16,N))*U.u;
LFP   = spm_int_J(P,M,U);
 
% and estimate spectral features under a MAR model
%--------------------------------------------------------------------------
mar = spm_mar(LFP,8);
mar = spm_mar_spectra(mar,M.Hz,1/U.dt);
CSD = abs(mar.P);
 
subplot(2,1,1)
plot([1:N]*U.dt,LFP)
xlabel('time')
title('LFP')
 
subplot(2,1,2)
plot(M.Hz,abs(CSD(:,1,1)),M.Hz,abs(CSD(:,1,2)),':')
xlabel('frequency')
title('[cross]-spectral density')
axis square
 
 
% inversion (in frequency space)
%==========================================================================
 
% data and confounds
%--------------------------------------------------------------------------
Y.y    = CSD;
nf     = size(Y.y,1);                     % number of frequency bins
Y.X0   = ones(nc*nc*nf,1);                % confounds (mean in log space)
Y.Q    = spm_Q(1/2,nf,1);                 % precision of noise AR(1/2)
 
% invert
%--------------------------------------------------------------------------
[Ep,Cp,S,F] = spm_nlsi_GN(M,[],Y);
 
 
 
% plot spectral density (after removing mean in log-space)
%==========================================================================
[G w] = spm_lfp_mtf(Ep,M);
 
% detrend
%------------------------------------------------------------------
g = spm_unvec(exp(spm_detrend(log(spm_vec(G)))),G);
y = spm_unvec(exp(spm_detrend(log(spm_vec(Y.y)))),G);
        
for i = 1:nc
    for j = 1:nc
        
        subplot(3,2,(i - 1)*nc + j)
        plot(w,g(:,i,j),w,y(:,i,j),':')
        title(sprintf('cross-spectral density %d,%d',i,j))
        xlabel('Power')
        axis square
        
        try axis(a),catch, a = axis; end
 
    end
end
legend({'predicted','observed'})
 
% plot parameters and estimates
%--------------------------------------------------------------------------
subplot(3,2,5)
bar(exp(spm_vec(P)))
title('true parameters')
 
subplot(3,2,6)
bar(exp(spm_vec(Ep)))
title('conditional expectation')
