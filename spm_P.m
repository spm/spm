function [P,p,Em,En,EN] = spm_P(c,k,Z,df,STAT,R,n,S)
% Returns the [un]corrected P value using unifed EC theory
% FORMAT [P p Em En EN] = spm_P(c,k,Z,df,STAT,R,n,S)
%
% c     - cluster number 
% k     - extent {RESELS}
% Z     - height {minimum over n values}
% df    - [df{interest} df{error}]
% STAT  - Statistical field
%		'Z' - Gaussian field
%		'T' - T - field
%		'X' - Chi squared field
%		'F' - F - field
%		'P' - Posterior probability
% R     - RESEL Count {defining search volume}
% n     - number of component SPMs in conjunction
% S     - Voxel count
%
% P     - corrected   P value  - P(n > kmax}
% p     - uncorrected P value  - P(n > k}
% Em    - expected total number of maxima {m}
% En    - expected total number of resels per cluster {n}
% EN    - expected total number of voxels {N}
%
%___________________________________________________________________________
%
% spm_P determines corrected and uncorrected p values, using the minimum
% of different valid methods. 
%
% See the individual methods for details
%
%     spm_P_RF
%     spm_P_Bonf
%
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Thomas Nichols
% $Id: spm_P.m 707 2006-12-06 16:42:20Z volkmar $


% set global var NOBONF to 1 to turn off Bonferroni
%---------------------------------------------------------------------------
global NOBONF; if ~isempty(NOBONF) & NOBONF, S = []; end

if (nargin < 8), S = []; end

[P,p,Em,En,EN] = spm_P_RF(c,k,Z,df,STAT,R,n);

% Use lower Bonferroni P value (if possible)
%===========================================================================
if ~isempty(S) & (c == 1 & k == 0) & ~(length(R) == 1 & R == 1)
    P = min(P,spm_P_Bonf(Z,df,STAT,S,n));
end
