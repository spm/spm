function [u] = spm_uc_Bonf(a,df,STAT,S,n)
% corrected critical height threshold at a specified significance level
% FORMAT [u] = spm_uc_Bonf(a,df,STAT,S,n)
% a     - critical probability - {alpha}
% df    - [df{interest} df{residuals}]
% STAT  - Statistical field
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
%       'P' - P - value
% S     - Voxel count
% n     - number of conjoint SPMs
%
% u     - critical height {corrected}
%
%___________________________________________________________________________
%
% spm_uc returns the corrected critical threshold at a specified significance
% level (a). If n > 1 a conjunction the probability over the n values of the 
% statistic is returned.
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Thomas Nichols
% $Id: spm_uc_Bonf.m 2690 2009-02-04 21:44:28Z guillaume $


u   = spm_u((a/S).^(1/n),df,STAT);
