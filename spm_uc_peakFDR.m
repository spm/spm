function [u, Ps] = spm_uc_peakFDR(q,df,STAT,R,n,Z,XYZ,ui)
% Peak False Discovery critical height threshold
% FORMAT [u] = spm_uc_peakFDR(q,df,STAT,R,n,Z,XYZ,ui)
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
% ui    - feature-inducing threshold
%
% u     - critical height threshold
% Ps    - Sorted p-values
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
% $Id: spm_uc_peakFDR.m 2975 2009-03-26 21:43:31Z guillaume $

ws       = warning('off','SPM:outOfRangePoisson');

% Extract list of local maxima whose height is above ui
%--------------------------------------------------------------------------
I        = find(Z >= ui);
Z        = Z(I);
XYZ      = XYZ(:,I);
[N, Z]   = spm_max(Z, XYZ);

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
