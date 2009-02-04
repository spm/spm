function [u] = spm_uc(a,df,STAT,R,n,S)
% corrected critical height threshold at a specified significance level
% FORMAT [u] = spm_uc(a,df,STAT,R,n,S)
% a     - critical probability - {alpha}
% df    - [df{interest} df{residuals}]
% STAT  - Statistical field
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
% R     - RESEL Count {defining search volume}
% n     - number of conjoint SPMs
% S     - Voxel count
%
% u     - critical height {corrected}
%
%__________________________________________________________________________
%
% spm_uc corrected critical thresholds, using the minimum of different
% valid methods.
%
% See the individual methods for details
%
%     spm_uc_RF
%     spm_uc_Bonf
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Thomas Nichols
% $Id: spm_uc.m 2690 2009-02-04 21:44:28Z guillaume $


% set global var NOBONF to 1 to turn off Bonferroni.
global NOBONF; if ~isempty(NOBONF) && NOBONF, S = []; end

if (nargin<6), S = []; end

u = spm_uc_RF(a,df,STAT,R,n);

if ~isempty(S)
  u = min(u,spm_uc_Bonf(a,df,STAT,S,n));
end
