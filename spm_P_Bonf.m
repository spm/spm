function P = spm_P_Bonf(Z,df,STAT,S,n)
% Return the corrected P value using Bonferroni
% FORMAT P = spm_P_Bonf(Z,df,STAT,S,n)
%
% Z     - height {minimum over n values}
% df    - [df{interest} df{error}]
% STAT  - Statistical field
%       'Z' - Gaussian field
%       'T' - T - field
%       'X' - Chi squared field
%       'F' - F - field
% S     - Voxel count
% n     - number of conjoint SPMs
%
% P     - corrected   P value  - P(STAT > Z)
%
%__________________________________________________________________________
%
% spm_P_Bonf returns the p-value of Z corrected by the Bonferroni
% inequality. 
%
% If n > 1 a conjunction probability over the n values of the statistic
% is returned.
%__________________________________________________________________________

% Thomas Nichols
% Copyright (C) 1999-2022 Wellcome Centre for Human Neuroimaging


P = S * spm_z2p(Z,df,STAT,n);
P = min(P,1);
