function [u] = spm_uc_RF(a,df,STAT,R,n)
% Corrected critical height threshold at a specified significance level
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
% level (a). If n > 1 a conjunction probability over the n values of the
% statistic is returned.
%__________________________________________________________________________
% Copyright (C) 1999-2012 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_uc_RF.m 4634 2012-02-01 19:01:09Z guillaume $


%-Find approximate value
%--------------------------------------------------------------------------
u  = spm_u((a/max(R))^(1/n),df,STAT);
du = 1e-6;

%-Approximate estimate using E{m}
%--------------------------------------------------------------------------
d  = 1;
while abs(d) > 1e-6
    [P, P, p] = spm_P_RF(1,0,u,df,STAT,R,n);
    [P, P, q] = spm_P_RF(1,0,u + du,df,STAT,R,n);
    d         = (a - p)/((q - p)/du);
    u         = u + d;
    if isinf(u), u=+Inf; return; end
end

%-Refined estimate using 1 - exp(-E{m})
%--------------------------------------------------------------------------
d  = 1;
while abs(d) > 1e-6
    p         = spm_P_RF(1,0,u,df,STAT,R,n);
    q         = spm_P_RF(1,0,u + du,df,STAT,R,n);
    d         = (a - p)/((q - p)/du);
    u         = u + d;
    if isinf(u), u=+Inf; return; end
end
