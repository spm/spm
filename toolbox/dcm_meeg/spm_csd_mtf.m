function [y,w,s,g] = spm_csd_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s,g] = spm_csd_mtf(P,M,U)
% FORMAT [y,w,s,g] = spm_csd_mtf(P,M)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% y - {y(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - directed transfer functions (complex)
% g - cross-spectral density (without channel noise)
%
% NB: requires M.u to specify the number of endogenous inputs
% This routine and will solve for the (hidden) steady state and use it as
% the expansion point for subsequent linear systems analysis (if trial
% specific effects are specified).
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_mtf.m 4807 2012-07-26 16:15:49Z guillaume $


% between-trial (experimental) inputs
%==========================================================================
try
    X = U.X;
    if ~size(X,1)
        X = sparse(1,0);
    end
catch
    
    % default inputs - one trial (no trial-specific effects)
    %----------------------------------------------------------------------
    X = sparse(1,0);
    
end

% compute log-spectral density
%==========================================================================

% frequencies of interest
%--------------------------------------------------------------------------
try
    dt = 1/(2*round(M.Hz(end)));
    N  = 1/dt;
    If = round(linspace(M.Hz(1),M.Hz(end),length(M.Hz)));
catch
    N  = 128;
    dt = 1/N;
    If = 1:N/2;
end
f    = (1:N/2)';                         % frequencies
w    = f(If);                            % frequencies selected

% number of channels and exogenous (neuronal) inputs or sources
%--------------------------------------------------------------------------
nc   = M.l;
ns   = length(M.u);

% spectrum of innovations (Gu) and noise (Gs and Gn)
%--------------------------------------------------------------------------
[Gu,Gs,Gn] = spm_csd_mtf_gu(P,M);


% cycle over trials (experimental conditions)
%==========================================================================
for  c = 1:size(X,1)
    
    % baseline parameters
    %----------------------------------------------------------------------
    Q  = P;
    
    % trial-specific effective connectivity
    %----------------------------------------------------------------------
    for i = 1:size(X,2)
        
        % extrinsic (forward and backwards) connections
        %------------------------------------------------------------------
        for j = 1:length(Q.A)
            Q.A{j} = Q.A{j} + X(c,i)*P.B{i};
        end
        
        % intrinsic connections
        %------------------------------------------------------------------
        try
            Q.H(:,1) = Q.H(:,1) + X(c,i)*diag(P.B{i});
        catch
            Q.G(:,1) = Q.G(:,1) + X(c,i)*diag(P.B{i});
        end
        
    end
    
    % solve for steady-state - for each condition
    %----------------------------------------------------------------------
    if nargin > 2
        M.x   = spm_dcm_neural_x(Q,M);
    end
    
    
    % delay operator - if parameterised
    %----------------------------------------------------------------------
    if nargout(M.f) == 3
        [unused,unused,D] = feval(M.f,M.x,M.u,Q,M);
    else
        D = 1;
    end
    
    
    % augment and bi-linearise
    %----------------------------------------------------------------------
    [M0,M1,L] = spm_bireduce(M,Q,D);
    
    % project onto spatial modes
    %----------------------------------------------------------------------
    if isfield(M,'U')
        L = M.U'*L;
    end
    
    % check for stability (a controllable system)
    %----------------------------------------------------------------------
    [E,V]  = eig(full(M0));
    V      = diag(V);
    if max(real(V)) > 0
        V  = real(V).*(real(V) < 0) + sqrt(-1)*imag(V);
        M0 = real(E*diag(V)*pinv(E));
    end
    
    % kernels
    %----------------------------------------------------------------------
    [K0,K] = spm_kernels(M0,M1,L,N,dt);
    
    % Transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    S      = fft(K);
    
    % [cross]-spectral density from neuronal innovations
    %----------------------------------------------------------------------
    G     = zeros(N/2,nc,nc);
    for i = 1:nc
        for j = 1:nc
            for k = 1:ns
                Gij      = S(:,i,k).*conj(S(:,j,k));
                Gij      = Gij((1:N/2) + 1).*Gu(:,k);
                G(:,i,j) = G(:,i,j) + Gij;
            end
        end
    end
    
    % save trial-specific frequencies of interest
    %----------------------------------------------------------------------
    g{c} = G(If,:,:);
    s{c} = S(If,:,:);
    
end

% and add channel noise
%--------------------------------------------------------------------------
for c = 1:length(g)
    
    G     = g{c};
    for i = 1:nc
        
        % channel specific noise
        %------------------------------------------------------------------
        G(:,i,i) = G(:,i,i) + Gs(If,i);
        
        % and cross-spectral density from common channel noise
        %------------------------------------------------------------------
        for j = 1:nc
            G(:,i,j) = G(:,i,j) + Gn(If);
        end
    end
    y{c} = G;
    
end

