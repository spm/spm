function [U] = spm_dcm_eeg_channelmodes(dipfit,Nm,xY)
% Returns the channel eigenmodes
% FORMAT [U] = spm_dcm_eeg_channelmodes(dipfit,Nm)
% FORMAT [U] = spm_dcm_eeg_channelmodes(dipfit,Nm,xY)
% dipfit  - spatial model specification
% Nm      - number of modes required (upper bound)
% xY      - data structure
% U       - channel eigenmodes
%__________________________________________________________________________
%
% Uses SVD (an eigensolution) to identify the patterns with the greatest 
% prior covariance; assuming independent source activity in the specified 
% spatial (forward) model. 
%
% if xY is specifed a CVA (a generalised eigensolution) will be used to 
% find the spatial modes that are best by the spatial model
%
% U is scaled to ensure trace(U'*L*L'*U) = Nm
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dcm_eeg_channelmodes.m 6806 2016-06-08 20:04:30Z karl $
 
% number of channels and modes
%--------------------------------------------------------------------------
if nargin < 2, Nm = 8; end

% Spatial modes
%--------------------------------------------------------------------------
[pE pC] = spm_L_priors(dipfit);

% evaluate eigenmodes of gain of covariance in sensor space
%--------------------------------------------------------------------------
dGdg  = spm_diff('spm_erp_L',pE,dipfit,1,'nocat');
L     = spm_cat(dGdg(find(spm_vec(pC))));

% eigen-mode reduction
%--------------------------------------------------------------------------
[U S] = spm_svd(L*L',exp(-8));
S     = diag(S);

if nargin > 2
    
    % reduce lead field to 32 or less d.f.
    %----------------------------------------------------------------------
    n    = min(size(U,2),32);
    L    = U(:,1:n)*diag(sqrt(S(1:n)));
    
    % response variable
    %----------------------------------------------------------------------
    for i = 1:length(xY.y)
        Y{i} = xY.y{i};
        Y{i} = Y{i}';
    end
    Y     = spm_cat(Y);
    CVA   = spm_cva(Y,L);
    U     = spm_en(CVA.w);
    S     = U'*(L*L')*U;
    S     = diag(S);
end

% eigen-mode reduction
%--------------------------------------------------------------------------
try
    U = U(:,1:Nm);
    S = S(  1:Nm);
end

% re-scale spatial projector
%--------------------------------------------------------------------------
U     = U/sqrt(mean(S));

