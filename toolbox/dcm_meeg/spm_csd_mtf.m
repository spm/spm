function [y,w,s] = spm_csd_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s] = spm_csd_mtf(P,M,U)
% FORMAT [y,w,s] = spm_csd_mtf(P,M)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects
%
% y - {y(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - directed transfer functions (complex)
%
% When called with U this function will return a cross-spectral respsone
% for each of the condition-specific paramters specifed in U.X; otherwise
% it simple return the complex CSD for the parameters in P.
%
% When the observer function M.g is pecifed the CSD repsonse is
% supplemented with channel noise in sneosr space; otherwsie the CSD
% pertins to hidden states.
%
% NB: requires M.u to specify the number of endogenous inputs
% This routine and will solve for the (hidden) steady state and use it as
% the expansion point for subsequent linear systems analysis (if trial
% specific effects are specified).
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_mtf.m 4814 2012-07-30 19:56:05Z karl $


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
    w  = round(linspace(M.Hz(1),M.Hz(end),length(M.Hz)));
catch
    N  = 128;
    dt = 1/N;
    w  = 1:N/2;
end

% number of channels and exogenous (neuronal) inputs or sources
%--------------------------------------------------------------------------
nc   = M.l;
ns   = length(M.u);
nw   = length(w);

% spectrum of innovations (Gu) and noise (Gs and Gn)
%--------------------------------------------------------------------------
if isfield(M,'g')
    [Gu,Gs,Gn] = spm_csd_mtf_gu(P,w);
else
    Gu         = spm_csd_mtf_gu(P,w);
end


% cycle over trials (experimental conditions)
%==========================================================================the
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
    
    
    % delay operator - if not specfied already and it is parameterised
    %----------------------------------------------------------------------
    if ~isfield(M,'D') && nargout(M.f) == 3
        [~,~,D] = feval(M.f,M.x,M.u,Q,M);
        M.D     = D;
    end
    
    
    % augment and bi-linearise
    %----------------------------------------------------------------------
    [M0,M1,L] = spm_bireduce(M,Q);

    % kernels
    %----------------------------------------------------------------------
    [~,K] = spm_kernels(M0,M1,L,N,dt);
    
    % Transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    S     = fft(K);
    S     = S(w + 1,:,:);
    
    % [cross]-spectral density from neuronal innovations
    %----------------------------------------------------------------------
    G     = zeros(nw,nc,nc);
    for i = 1:nc
        for j = 1:nc
            for k = 1:ns
                Gij      = S(:,i,k).*conj(S(:,j,k));
                G(:,i,j) = G(:,i,j) + Gij.*Gu(:,k);
            end
        end
    end
    
    % save trial-specific frequencies of interest
    %----------------------------------------------------------------------
    g{c}  = G;
    s{c}  = S;
    
end

% and add channel noise
%==========================================================================
if isfield(M,'g')
    
    for c = 1:length(g)
        G = g{c};
        for i = 1:nc
            
            % channel specific noise
            %--------------------------------------------------------------
            try
                G(:,i,i) = G(:,i,i) + Gs(:,i);
            catch
                G(:,i,i) = G(:,i,i) + Gs(:,1);
            end
            
            % and cross-spectral density from common channel noise
            %--------------------------------------------------------------
            for j = 1:nc
                G(:,i,j) = G(:,i,j) + Gn;
            end
        end
        y{c} = G;
        
    end
else
    y = g;
end

