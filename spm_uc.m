function [u] = spm_uc(a,df,STAT,R,n,S)
% corrected critical height threshold at a specified significance level
% FORMAT [u] = spm_uc(a,df,STAT,R,n,S)
% a     - critical probability - {alpha}
% df    - [df{interest} df{residuals}]
% STAT  - Statistical field
%		'Z' - Gaussian field
%		'T' - T - field
%		'X' - Chi squared field
%		'F' - F - field
% R     - RESEL Count {defining search volume}
% n     - number of conjoint SPMs
% S     - Voxel count
%
% u     - critical height {corrected}
%
%___________________________________________________________________________
%
% spm_uc corrected critical thresholds, using the minimum of different
% valid methods.
%
% See the individual methods for details
%
%     spm_uc_RF
%     spm_uc_Bonf
%
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Thomas Nichols
% $Id: spm_uc.m 707 2006-12-06 16:42:20Z volkmar $


% set global var NOBONF to 1 to turn off Bonferroni.
global NOBONF; if ~isempty(NOBONF) & NOBONF, S = []; end

if (nargin<6), S = []; end

u = spm_uc_RF(a,df,STAT,R,n);

if ~isempty(S)
  u = min(u,spm_uc_Bonf(a,df,STAT,S,n));
end
