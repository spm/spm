function [y,H] = spm_bsd_int_parfor(P,M,U)
% Spectral response of a NMM (transfer function x noise spectrum)
% FORMAT [y,w,s,g] = spm_csd_mtf(P,M,U)
% FORMAT [y,w,s,g] = spm_csd_mtf(P,M)
%
% P - parameters
% M - neural mass model structure
% U - trial-specific effects (induces expansion around steady state)
%
% y - {y(N,nc,nc}} - cross-spectral densityDCM for nc channels {trials}
%                  - for N frequencies in M.Hz [default 1:64Hz]
% w - frequencies
% s - modulation transfer functions (complex)
% g - normalised modulation transfer function (true Granger causality)
%
% When called with U this function will return a cross-spectral response
% for each of the condition-specific parameters specified in U.X; otherwise
% it returns the complex CSD for the parameters in P (using the expansion
% point supplied in M.x)
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
%__________________________________________________________________________

% Johan Medrano
% Copyright (C) 2005-2022 Wellcome Centre for Human Neuroimaging


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
if isfield(M,'Hz')
    w    = M.Hz;
else
    w    = 1:64;
    M.Hz = w;
end

L  = spm_lx_erp(P,M.dipfit);

% project onto spatial modes
%--------------------------------------------------------------------------
if isfield(M,'U')
    L = M.U'*L;
end

if nargout > 1
    H = {}; 
end

% cycle over trials (experimental conditions)
%==========================================================================
nargs = nargout; 
parfor  c = 1:size(X,1)
    

    % condition-specific parameters
    %----------------------------------------------------------------------
%     Q   = spm_gen_Q(P,X(c,:));
%     Q = spm_unvec(spm_vec(P)
    if numel(P.B)
        Qp = keepfields(P, fieldnames(P.B{1})); 
        Q = spm_vec(Qp);
        for i = 1:size(X, 2)
            Q = Q + X(c,i)*spm_vec(P.B{i});
        end
        Q = spm_unvec(Q, Qp); 
    else
        Q = P; 
    end
    if nargs > 1
        [s, Hc] = spm_bsd_csd(Q, M); 
        for j = 1:numel(Hc.modes)
            Hc.modes{j} = source2sensor(Hc.modes{j}, L); 
        end
        Hc.noise = source2sensor(Hc.noise, L); 
        H{c}   = Hc; 
    else 
        s  = spm_bsd_csd(Q, M); 

    end
    y{c} = source2sensor(s, L);     
end
end

function y = source2sensor(s, L)
    % if nargin < 3
    %     norm = sum(s, 1); 
    % end  
    % s = s/norm;
    norm = sqrt(sum(L.^2, 1));
    if sqrt(sum(L.^2, 1)) > 0
        L = L./norm;
    end
    y = zeros(size(s, 1), size(L, 1), size(L,1)); 
    for ii = 1:size(s, 1)
        si = reshape(s(ii, :, :), size(s,2,3));
        y(ii, :, :) = L * si * L'; 
    end
end

