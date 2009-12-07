function [u] = spm_uc_RF(a,df,STAT,R,n)
% corrected critical height threshold at a specified significance level
% FORMAT [u] = spm_uc_RF(a,df,STAT,R,n)
% a     - critical probability - {alpha}
% df    - [df{interest} df{residuals}]
% STAT  - Statistical field
%         'Z' - Gaussian field
%         'T' - T field
%         'X' - Chi-squared field
%         'F' - F field
% R     - RESEL Count {defining search volume}
% n     - number of conjoint SPMs
%
% u     - critical height {corrected}
%
%__________________________________________________________________________
%
% spm_uc returns the corrected critical threshold at a specified significance
% level (a). If n > 1 a conjunction the probability over the n values of
% the statistic is returned.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_uc_RF.m 3614 2009-12-07 18:04:24Z guillaume $


% find approximate value
%--------------------------------------------------------------------------
u  = spm_u((a/max(R))^(1/n),df,STAT);
du = 1e-6;

% approximate estimate using E{m}
%--------------------------------------------------------------------------
d  = 1;
while abs(d) > 1e-6
    [P P p] = spm_P_RF(1,0,u,df,STAT,R,n);
    [P P q] = spm_P_RF(1,0,u + du,df,STAT,R,n);
    d       = (a - p)/((q - p)/du);
    u       = u + d;
end

% refined estimate using 1 - exp(-E{m})
%--------------------------------------------------------------------------
d  = 1;
while abs(d) > 1e-6
    p       = spm_P_RF(1,0,u,df,STAT,R,n);
    q       = spm_P_RF(1,0,u + du,df,STAT,R,n);
    d       = (a - p)/((q - p)/du);
    u       = u + d;
end
