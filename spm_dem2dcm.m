function [DCM]   = spm_dem2dcm(DEM,DCM)
% reorganisation of posteriors and priors into DCM format
% FORMAT [DCM]   = spm_dem2dcm(DEM)
% FORMAT [DEM]   = spm_dem2dcm(DEM,DCM)
%
% DEM - structure array (hierarchicial model)
% DCM - structure array (flat model)
%
% -------------------------------------------------------------------------
%     DCM.M.pE - prior expectation of parameters
%     DCM.M.pC - prior covariances of parameters
%     DCM.Ep   - posterior expectations
%     DCM.Cp   - posterior covariance
%     DCM.F   - free energy
%
% For hierarchical models (DEM.M) the first level with non-zero prorior
% varaince on the paramters will be extracted.
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dem2dcm.m 6508 2015-07-25 15:23:25Z karl $


% check for arrays
%--------------------------------------------------------------------------
if iscell(DEM)
    for i = 1:size(DEM,1)
        for j = 1:size(DEM,2)
            if nargin < 2
                DCM{i,j}   = spm_dem2dcm(DEM{i,j});
            else
                DCM{i,j}   = spm_dem2dcm(DEM{i,j},DCM{i,j});
            end
        end
    end
    return
end


% get level
%--------------------------------------------------------------------------
for j = 1:length(DEM.M)
    if any(any(DEM.M(j).pC)), break, end
end

% get indices for covariance
%--------------------------------------------------------------------------
k     = spm_length(DEM.qP.P(j));
k     = spm_length(DEM.qP.P(1:(j - 1))) + (1:k);

% re-organise
%--------------------------------------------------------------------------
if nargin < 2
    
    DCM.M.pE      = DEM.M(j).pE;
    DCM.M.pC      = DEM.M(j).pC;
    DCM.Ep        = DEM.qP.P{j};
    DCM.Cp        = DEM.qP.C(k,k);
    DCM.F         = DEM.F(end);
    
else
    
    DEM.M(j).pE   = spm_unvec(DCM.M.pE,DEM.M(j).pE);
    DEM.M(j).pC   = DCM.M.pC;
    DEM.qP.P{j}   = spm_unvec(DCM.Ep,DEM.qP.P{j});
    DEM.qP.C(k,k) = DCM.Cp;
    DEM.F         = DCM.F;
    
    % output argument
    %----------------------------------------------------------------------
    DCM = DEM;
end





