function [u] = spm_u(a,df,STAT)
% uncorrected critical height threshold at a specified significance level
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
% significance
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_u.m 2690 2009-02-04 21:44:28Z guillaume $


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
