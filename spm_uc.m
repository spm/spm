function [u] = spm_uc(a,df,STAT,R,n,S)
% corrected critical height threshold at a specified significance level
% FORMAT [u] = spm_uc(a,df,STAT,R,n,S)
% a     - critical probability - {alpha}
% df    - [df{interest} df{residuals}]
% STAT  - Statisical feild
%		'Z' - Gaussian feild
%		'T' - T - feild
%		'X' - Chi squared feild
%		'F' - F - feild
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
% $Id$


% set global var NOBONF to 1 to turn off Bonferroni.
global NOBONF; if ~isempty(NOBONF) & NOBONF, S = []; end

if (nargin<6), S = []; end

u = spm_uc_RF(a,df,STAT,R,n);

if ~isempty(S)
  u = min(u,spm_uc_Bonf(a,df,STAT,S,n));
end
