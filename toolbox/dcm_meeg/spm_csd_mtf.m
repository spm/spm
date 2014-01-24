function [y,w,s] = spm_csd_mtf(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s] = spm_csd_mtf(P,M,U)
% FORMAT [y,w,s] = spm_csd_mtf(P,M)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects (induces expansion around steady state)
%
% y - {y(N,nc,nc}} - cross-spectral density for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - modulation transfer functions (complex)
%
% When called with U this function will return a cross-spectral response
% for each of the condition-specific parameters specified in U.X; otherwise
% it returns the complex CSD for the parameters in P.
%
% When the observer function M.g is specified the CSD response is
% supplemented with channel noise in sensor space; otherwise the CSD
% pertains to hidden states.
%
% NB: requires M.u to specify the number of endogenous inputs
% This routine and will solve for the (hidden) steady state and use it as
% the expansion point for subsequent linear systems analysis (if trial
% specific effects are specified).
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_csd_mtf.m 5853 2014-01-24 20:38:11Z karl $



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
    w    = round(linspace(M.Hz(1),M.Hz(end),length(M.Hz)));
catch
    w    = 1:64;
    M.Hz = w;
end

% number of channels and exogenous (neuronal) inputs or sources
%--------------------------------------------------------------------------
nc   = M.l;
nw   = length(M.Hz);

% spectrum of innovations (Gu) and noise (Gs and Gn)
%--------------------------------------------------------------------------
if isfield(M,'g')
    [Gu,Gs,Gn] = spm_csd_mtf_gu(P,w);
else
    Gu         = spm_csd_mtf_gu(P,w);
end


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
            
            % CMM-NMDA specific modulation on Extrinsic NMDA connections
            %--------------------------------------------------------------
            try
                Q.AN{j} = Q.AN{j} + X(c,i)*P.BN{i};
            end
        end
        
        % intrinsic connections
        %------------------------------------------------------------------
        Q.G(:,1) = Q.G(:,1) + X(c,i)*diag(P.B{i});
        
    end
    
    % solve for steady-state - if exogenous inputs are specified
    %----------------------------------------------------------------------
    if nargin > 2
        M.x = spm_dcm_neural_x(Q,M);
    end
    
    % transfer functions (FFT of kernel)
    %----------------------------------------------------------------------
    S     = spm_dcm_mtf(Q,M);
    

    % predicted cross-spectral density
    %----------------------------------------------------------------------
    G     = zeros(nw,nc,nc);
    for i = 1:nw
        G(i,:,:) = squeeze(S(i,:,:))*diag(Gu(i,:))*squeeze(S(i,:,:))';
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
