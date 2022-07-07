function [u] = spm_u(a,df,STAT)
% Uncorrected critical height threshold at a specified significance level
% FORMAT [u] = spm_u(a,df,STAT)
% a     - critical probability - {alpha}
% df    - [df{interest} df{error}]
% STAT  - Statistical field
%               'Z' - Gaussian field
%               'T' - T - field
%               'X' - Chi squared field
%               'F' - F - field
%               'P' - P - value
%
% u     - critical height {uncorrected}
%__________________________________________________________________________
%
% spm_u returns the uncorrected critical threshold at a specified 
% significance.
%__________________________________________________________________________

% Karl Friston
% Copyright (C) 1995-2022 Wellcome Centre for Human Neuroimaging


if     STAT == 'Z'

    u   = spm_invNcdf(1 - a      );

elseif STAT == 'T'

    u   = spm_invTcdf(1 - a,df(2));

elseif STAT == 'X'

    u   = spm_invXcdf(1 - a,df(2));

elseif STAT == 'F'

    u   = spm_invFcdf(1 - a,df   );

elseif STAT == 'P'

    u   = a;

end
