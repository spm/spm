function [u, Ps] = spm_uc_peakFDR(q,df,STAT,R,n,Z,XYZ,ui,G)
% Peak False Discovery critical height threshold
% FORMAT [u, Ps] = spm_uc_peakFDR(q,df,STAT,R,n,Z,XYZ,ui[,G])
%
% q     - prespecified upper bound on False Discovery Rate
% df    - [df{interest} df{residuals}]
% STAT  - statistical field
%         'Z' - Gaussian field
%         'T' - T - field
%         'X' - Chi squared field
%         'F' - F - field
% R     - RESEL Count {defining search volume}
% n     - conjunction number
% Z     - height {minimum over n values}
%         or mapped statistic image(s)
% XYZ   - locations [x y x]' {in voxels}
%         or vector of indices of elements within mask
%         or mapped mask image
% ui    - feature-inducing threshold
% G     - patch structure (for surface-based inference)
%
% u     - critical height threshold
% Ps    - sorted p-values
%__________________________________________________________________________
%
% References
%
% J.R. Chumbley and K.J. Friston, "False discovery rate revisited: FDR and 
% topological inference using Gaussian random fields". NeuroImage,
% 44(1):62-70, 2009.
%
% J.R. Chumbley, K.J. Worsley, G. Flandin and K.J. Friston, "Topological
% FDR for NeuroImaging". NeuroImage, 49(4):3057-3064, 2010.
%__________________________________________________________________________
% Copyright (C) 2009-2013 Wellcome Trust Centre for Neuroimaging

% Justin Chumbley & Guillaume Flandin
% $Id: spm_uc_peakFDR.m 5224 2013-02-01 12:19:17Z guillaume $


ws = warning('off','SPM:outOfRangePoisson');

% Read statistical value from disk if needed
%--------------------------------------------------------------------------
if isstruct(Z)
    Vs         = Z;
    Vm         = XYZ;
    [Z, XYZmm] = spm_read_vols(Vs(1),true);
    for i=2:numel(Vs)
        Z      = min(Z, spm_read_vols(Vs(i),true));
    end
    Z          = Z(:)';
    XYZ        = round(Vs(1).mat \ [XYZmm; ones(1, size(XYZmm, 2))]);
    try
        M      = spm_vol(spm_file(Vs(1).fname,'basename','mask'));
        I      = spm_get_data(M,XYZ,false) > 0;
    catch
        I      = ~isnan(Z) & Z~=0;
    end
    XYZ        = XYZ(1:3,I);
    Z          = Z(I);
    if isstruct(Vm)
        Vm     = spm_get_data(Vm,XYZ,false) > 0;
    end

    XYZ        = XYZ(:,Vm);
    Z          = Z(:,Vm);
end

% Extract list of local maxima whose height is above ui
%--------------------------------------------------------------------------
I        = find(Z >= ui);
Z        = Z(I);
XYZ      = XYZ(:,I);
if nargin == 8
    [N,Z] = spm_max(Z, XYZ);
else
    [N,Z] = spm_mesh_max(Z, XYZ, G);
end

% Expected Euler characteristic for level ui
%--------------------------------------------------------------------------
[P,p,Eu] = spm_P_RF(1,0,ui,df,STAT,R,n);

% Expected Euler characteristic for level Z(i)
%--------------------------------------------------------------------------
Ez       = zeros(1,numel(Z));
for i = 1:length(Z)
    [P,p,Ez(i)] = spm_P_RF(1,0,Z(i),df,STAT,R,n);
end

% Uncorrected p-value for peaks using Random Field Theory
%--------------------------------------------------------------------------
[Ps, J]  = sort(Ez ./ Eu, 'ascend');

S        = length(Ps);

% Calculate FDR inequality RHS
%--------------------------------------------------------------------------
cV       = 1;    % Benjamini & Yeuketeli cV for independence/PosRegDep case
Fi       = (1:S)/S*q/cV;

% Find threshold
%--------------------------------------------------------------------------
I        = find(Ps <= Fi, 1, 'last');
if isempty(I)
    u    = Inf;
else
    u    = Z(J(I));
end

warning(ws);
