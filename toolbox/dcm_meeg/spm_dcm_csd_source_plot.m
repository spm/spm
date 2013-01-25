function [G] = spm_dcm_csd_source_plot(model,s,P,N)
% Spectral response (G) of a single source neural mass model
% FORMAT [G] = spm_dcm_csd_source_plot(model,s)
%
% model - 'ERP', 'SEP', 'CMC', 'LFP', 'NMM' or 'MFM'
% s     - indices of hidden neuronal states to plot
% P     - parameters
% N     - twice the maximum frequency
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_csd_source_plot.m 5210 2013-01-25 15:31:46Z guillaume $


% Create model
%==========================================================================
try, s; catch, s = [1 2 3 7]; end
try, N; catch, N = 256;       end


% intial states and equations of motion
%--------------------------------------------------------------------------
[x,f] = spm_dcm_x_neural(P,model);

% create DCM
%--------------------------------------------------------------------------
M.f   = f;
M.x   = x;
M.n   = length(spm_vec(x));
M.m   = size(P.C,1);

% solve for steady state
%--------------------------------------------------------------------------
M.x   = spm_dcm_neural_x(P,M);

% exogenous (neuronal) inputs
%--------------------------------------------------------------------------
M.u   = 0;


% compute spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
dt   = 1/N;
f    = (1:N/2)';

% spectrum of innovations (Gu)
%--------------------------------------------------------------------------
Gu   =  f.^(-1/4);

% get delay operator, augment and bi-linearise (with delays)
%--------------------------------------------------------------------------
[M0,M1,D] = feval(M.f,M.x,M.u,P,M);
M.D       = D;
[M0,M1,L] = spm_bireduce(M,P);

% compute modulation transfer function using FFT of the kernels
%--------------------------------------------------------------------------
[K0,K1]   = spm_kernels(M0,M1,L,N,dt);
[N,nc,nu] = size(K1);


% [cross]-spectral density
%--------------------------------------------------------------------------
ns    = length(s);
G     = zeros(N/2,ns,ns);
for i = 1:ns
    for j = i:ns
        for k = 1:nu
            Si       = fft(K1(:,s(i),k));
            Sj       = fft(K1(:,s(j),k));
            Gij      = Si.*conj(Sj);
            Gij      = Gij((1:N/2) + 1).*Gu;
            G(:,i,j) = G(:,i,j) + Gij;
        end
        
        % Graphics
        %==================================================================
        subplot(ns,ns,(i - 1)*ns + j)
        plot(f,abs(G(:,i,j)))
        xlabel('frequency')
        axis square
        
    end
end

drawnow
