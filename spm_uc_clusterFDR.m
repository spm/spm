function [u, Ps, ue] = spm_uc_clusterFDR(q,df,STAT,R,n,Z,XYZ,V2R,ui)
% Cluster False Discovery critical extent threshold
% FORMAT [u, Ps, ue] = spm_uc_clusterFDR(q,df,STAT,R,n,Z,XYZ,ui)
%
% q     - Prespecified upper bound on False Discovery Rate
% df    - [df{interest} df{residuals}]
% STAT  - Statistical field
%         'Z' - Gaussian field
%         'T' - T - field
%         'X' - Chi squared field
%         'F' - F - field
% R     - RESEL Count {defining search volume}
% n     - Conjunction number
% Z     - height {minimum over n values}
% XYZ   - locations [x y x]' {in voxels}
% V2R   - voxel to resel
% ui    - feature-inducing threshold
%
% u     - critical extent threshold
% Ps    - Sorted p-values
% ue    - critical extent threshold for FWE
%__________________________________________________________________________
%
% References
%
% J.R. Chumbley and K.J. Friston, "False discovery rate revisited: FDR and 
% topological inference using Gaussian random fields". NeuroImage,
% 44(1):62-70, 2009.
%
% J.R. Chumbley, K.J. Worsley, G. Flandin and K.J. Friston, "Topological
% FDR for NeuroImaging". Under revision.
%__________________________________________________________________________
% Copyright (C) 2009 Wellcome Trust Centre for Neuroimaging

% Justin Chumbley & Guillaume Flandin
% $Id: spm_uc_clusterFDR.m 2975 2009-03-26 21:43:31Z guillaume $

% Threshold the statistical field 
%--------------------------------------------------------------------------
XYZ      = XYZ(:,Z >= ui);

% Extract size of excursion sets 
%--------------------------------------------------------------------------
N        = spm_clusters(XYZ);
if ~isempty(N)
    N    = histc(N,(0:max(N))+0.5);
end
N        = N(1:end-1);
N        = N .* V2R;

% Compute uncorrected p-values based on N using Random Field Theory
%--------------------------------------------------------------------------
Ps       = zeros(1,numel(N));
Pk       = zeros(1,numel(N));
ws       = warning('off','SPM:outOfRangePoisson');
for i = 1:length(N)
    [Pk(i), Ps(i)] = spm_P_RF(1,N(i),ui,df,STAT,R,n);
end
warning(ws);
[Ps, J]  = sort(Ps, 'ascend');

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
    u    = N(J(I)) / V2R;
end

% As we're there, also determine the FWE critical extent threshold
%--------------------------------------------------------------------------
if nargout == 3
    [Pk, J]  = sort(Pk, 'ascend');
    I        = find(Pk <= q, 1, 'last');
    if isempty(I)
        ue   = Inf;
    else
        ue   = N(J(I)) / V2R;
    end
end