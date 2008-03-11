function [y,w] = spm_lfp_mtf_laplace(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [G,w] = spm_lfp_mtf(P,M,U)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% G - {G(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_lfp_mtf_laplace.m 1199 2008-03-11 15:11:02Z rosalyn $


% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
try
    dt   = 1/(2*round(M.Hz(end)));
    N    = 1/dt;
    If   = round(linspace(M.Hz(1),M.Hz(end),length(M.Hz)));
catch
    N  = 128;
    dt = 1/N;
    If = 1:N/2;
end
f    = [1:N/2]';
w    = M.Hz*2*pi;


% spectrum of innovations (Gu)
%--------------------------------------------------------------------------
Gu   = exp(P.a)*M.Hz.^(-1)*2;             % spectral density of (AR) input
Gu   = Gu + exp(P.b);                     % spectral density of IID input

% channel noise (per channel)
%--------------------------------------------------------------------------

Gn   = exp(P.c)/8;

% trial-specific effects
%==========================================================================
try, X = U.X; catch, X = sparse(1,0); end

% basline parameters
%--------------------------------------------------------------------------
Q  = P;

% cycle over trials
%--------------------------------------------------------------------------
for  c = 1:size(X,1)
    % exogenous (neuronal) inputs
    %----------------------------------------------------------------------
    M.u  = sparse(length(P.C),1);

    % trial-specific effective connectivity
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        P.A{1} = Q.A{1} + X(c,i)*P.B{i};        % forward   connections
        P.A{2} = Q.A{2} + X(c,i)*P.B{i};        % backward  connections
        P.A{3} = Q.A{3} + X(c,i)*P.B{i};        % lateral   connections

        P.H    = Q.H + X(c,i)*diag(P.B{i});      % intrinsic connections
    end

    % augment and bi-linearise
    %----------------------------------------------------------------------
    [M0,M1,L] = spm_bireduce(M,P);
    % project onto spatial modes
    %----------------------------------------------------------------------
    try
        L = M.U'*L;
    end

    % compute modulation transfer function using FFT of the kernels
    %----------------------------------------------------------------------
    [K0,K1]   = spm_kernels(M0,M1,L,N,dt);


    % [cross]-spectral density
    %--------------------------------------------------------------------------
    [N,nc,nu] = size(K1);
    G     = zeros(length(w),nc,nc);

    A = full(M0(2:end,2:end));
    for i = 1:nc
        B(:,i) = full(M1{i}(2:end,1));
    end

    C = full(L(:,2:end));
    D = zeros(nc,nu);

    for i = 1:nc
        for j = 1:nc
            % exogenous components
            %-----------------------------------------------------------
            for k =1:nu
                [bi,ai]     =  spm_ss2tf(A,B,C(i,:),D(i,:),k);
                [bj,aj]     =  spm_ss2tf(A,B,C(j,:),D(j,:),k);
                Si          =  spm_freqs(bi,ai,w);
                Sj          =  spm_freqs(bj,aj,w);
                Gij         =  Si.*conj(Sj);
                Gij         = (abs(Gij).*Gu)';

                G(:,i,j) = G(:,i,j) + Gij;
                G(:,j,i) = G(:,j,i) + Gij;

            end
            % endogenous components
            %------------------------------------------------------------
            if i == j
                G(:,j,i)=G(:,j,i)+Gn;
            end
                
        end
    end

    % save frequencies of interest
    %----------------------------------------------------------------------
    y{c} = G;

end
