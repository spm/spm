function [DCM]   = spm_dem2dcm(DEM,j)
% reorganisation of posterior is and priors into DCM format
% FORMAT [DCM]   = spm_dem2dcm(DEM,[j])
% FORMAT [DEM]   = spm_dem2dcm(DEM,[j],DCM)
%
% DEM - structure array
% DCM - structure array 
% j   - level to extract/replace [all levels]
%
% -------------------------------------------------------------------------
%     DCM.M.pE - prior expectation of parameters
%     DCM.M.pC - prior covariances of parameters
%     DCM.Ep   - posterior expectations
%     DCM.Cp   - posterior covariance
%     DCM.F   - free energy
%__________________________________________________________________________
% Copyright (C) 2005 Wellcome Trust Centre for Neuroimaging
 
% Karl Friston
% $Id: spm_dem2dcm.m 6506 2015-07-24 10:26:51Z karl $

% check input arguments
%--------------------------------------------------------------------------
if nargin < 2, j = 1:length(DEM.M); end
if isempty(j), j = 1:length(DEM.M); end
if nargin < 3
    
    
    
    DEM   = P;
    j     = 2;
    k     = spm_length(DEM{1}.qP.P(j));
    k     = spm_length(DEM{1}.qP.P(1:(j - 1))) + (1:k);
    
    DCM.M.pE = DEM.M(j).pE;
    DCM.M.pC = DEM.M(j).pC;
    DCM.Ep   = DEM.qP.P{j};
    DCM.Cp   = DEM.qP.C(k,k);
    DCM.F    = DEM.F(end);
else
    
    
    j     = 2;
    k     = spm_length(DEM{1}.qP.P(j));
    k     = spm_length(DEM{1}.qP.P(1:(j - 1))) + (1:k);
    for i = 1:numel(P)
        DEM{i}.M(j).pE   = P{i}.M.pE;
        DEM{i}.M(j).pC   = P{i}.M.pC;
        DEM{i}.qP.P{j}   = P{i}.Ep;
        DEM{i}.qP.C(k,k) = P{i}.Cp;
        DEM{i}.F         = P{i}.F;
    end
    P     = DEM;
    
    
end





