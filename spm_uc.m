function [u] = spm_uc(a,df,STAT,R,n,S)
% Corrected critical height threshold at a specified significance level
% FORMAT [u] = spm_uc(a,df,STAT,R,n,S)
% a     - critical probability - {alpha}
% df    - [df{interest} df{residuals}]
% STAT  - Statistical field
%         'Z' - Gaussian field
%         'T' - T - field
%         'X' - Chi squared field
%         'F' - F - field
% R     - RESEL Count {defining search volume}
% n     - number of conjoint SPMs
% S     - Voxel count
%
% u     - critical height {corrected}
%__________________________________________________________________________
%
% spm_uc returns the corrected critical height threshold at a specified
% significance level (a), using the minimum of different valid methods.
%
% See also: spm_uc_RF, spm_uc_Bonf
%__________________________________________________________________________

% Thomas Nichols
% Copyright (C) 2001-2022 Wellcome Centre for Human Neuroimaging


if nargin < 6, S = []; end

u = spm_uc_RF(a,df,STAT,R,n);

% Compare with Bonferroni corrected height threshold (if possible)
%--------------------------------------------------------------------------
if ~isempty(S) && ~isequal(R, 1)
  u = min(u, spm_uc_Bonf(a,df,STAT,S,n));
end
