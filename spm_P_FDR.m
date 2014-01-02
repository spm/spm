function [P] = spm_P_FDR(Z,df,STAT,n,Ps)
% Returns the corrected FDR P value
%
% FORMAT [P] = spm_P_FDR(Z,df,STAT,n,Ps)
%
% Z     - height (minimum of n statistics)
% df    - [df{interest} df{error}]
% STAT  - Statistical field
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
%       'P' - P - value
% n     - Conjunction number
% Ps    - Vector of sorted (ascending) p-values in search volume
%
% P     - corrected FDR P value
%
%
% FORMAT [P] = spm_P_FDR(p)
%
% p     - vector or array of all uncorrected P values, from which
%         non-finite values will be excluded (but zeros and ones are kept)
%
% P     - corrected FDR P values (a vector or array of the same shape as p)
%
%__________________________________________________________________________
%
% The Benjamini & Hochberch (1995) False Discovery Rate (FDR) procedure
% finds a threshold u such that the expected FDR is at most q.  spm_P_FDR
% returns the smallest q such that Z>u. 
%
% Background
%
% For a given threshold on a statistic image, the False Discovery Rate
% is the proportion of suprathreshold voxels which are false positives.
% Recall that the thresholding of each voxel consists of a hypothesis
% test, where the null hypothesis is rejected if the statistic is larger
% than threshold.  In this terminology, the FDR is the proportion of
% rejected tests where the null hypothesis is actually true.
%
% A FDR proceedure produces a threshold that controls the expected FDR
% at or below q.  The FDR adjusted p-value for a voxel is the smallest q
% such that the voxel would be suprathreshold.
%
% In comparison, a traditional multiple comparisons proceedure
% (e.g. Bonferroni or random field methods) controls Familywise Error
% rate (FWER) at or below alpha.  FWER is the *chance* of one or more
% false positives anywhere (whereas FDR is a *proportion* of false
% positives).  A FWER adjusted p-value for a voxel is the smallest alpha
% such that the voxel would be suprathreshold. 
%
% If there is truely no signal in the image anywhere, then a FDR
% proceedure controls FWER, just as Bonferroni and random field methods
% do. (Precisely, controlling E(FDR) yields weak control of FWE).  If
% there is some signal in the image, a FDR method should be more powerful
% than a traditional method.
%
% For careful definition of FDR-adjusted p-values (and distinction between
% corrected and adjusted p-values) see Yekutieli & Benjamini (1999).
%
%
% References
%
% Benjamini & Hochberg (1995), "Controlling the False Discovery Rate: A
% Practical and Powerful Approach to Multiple Testing". J Royal Stat Soc,
% Ser B.  57:289-300.
%
% Benjamini & Yekutieli (2001), "The Control of the false discovery rate
% in multiple testing under dependency". Annals of Statistics, 
% 29(4):1165-1188.
%
% Yekutieli & Benjamini (1999). "Resampling-based false discovery rate
% controlling multiple test procedures for correlated test
% statistics".  J of Statistical Planning and Inference, 82:171-196.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Thomas Nichols and Ged Ridgway
% $Id: spm_P_FDR.m 5824 2014-01-02 14:50:13Z guillaume $

%-Set Benjamini & Yeuketeli cV for independence/PosRegDep case
%--------------------------------------------------------------------------
cV = 1; 

if nargin == 1
    %-FORMAT [P] = spm_P_FDR(p), with p -> Z, (all uncorrected p-values)
    %======================================================================
    
    %-Exclude non-finite (allowing use on a masked image) and sort
    %----------------------------------------------------------------------
    mask    = isfinite(Z);
    Z       = Z(mask);
    S       = numel(Z);
    [Z, si] = sort(Z(:));
    [i, ui] = sort(si); % "unsorting indices"
    
    %-"Corrected" p-values
    %----------------------------------------------------------------------
    Z       = Z*S./(1:S)'*cV;
    
    %-"Adjusted" p-values
    %----------------------------------------------------------------------
    Z(end + 1) = 1; % (sentinel 1 at end, dropped in ui indexing below)
    for i = S:-1:1
        Z(i) = min(Z([i i+1]));
    end
    
    %-Unsort, unmask and return
    %----------------------------------------------------------------------
    Z = Z(ui);
    P = nan(size(mask));
    P(mask) = Z;
    return
end


%-FORMAT [P] = spm_P_FDR(Z,df,STAT,n,Ps)
%==========================================================================

%-Calculate p-value of Z
%--------------------------------------------------------------------------
PZ   = spm_z2p(Z,df,STAT,n);

%-Calculate FDR p-values
%--------------------------------------------------------------------------
% If Z is a value in the statistic image, then the adjusted p-value
% defined in Yekutieli & Benjamini (1999) (eqn 3) is obtained.  If Z
% isn't a value in the image, then the adjusted p-value for the next
% smallest statistic value (next largest uncorrected p) is returned.

%-"Corrected" p-values
%--------------------------------------------------------------------------
S     = length(Ps);
Qs    = Ps*S./(1:S)'*cV;    

%-"Adjusted" p-values
%--------------------------------------------------------------------------
P     = zeros(size(PZ));
for i = 1:numel(PZ)

    %-Find PZ(i) in Ps, or smallest Ps(j) s.t. Ps(j) >= PZ(i)
    %----------------------------------------------------------------------
    I = find(Ps>=PZ(i), 1 );  

    %-"Adjusted" p-values
    %----------------------------------------------------------------------
    if isempty(I)
        P(i) = 1;
    else
        P(i) = min(Qs(I:S));
    end
end
