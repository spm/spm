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
%
% u     - critical height {uncorrected}
%___________________________________________________________________________
% spm_u returns the uncorrected critical threshold at a specified significance
%
%___________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_u.m 1143 2008-02-07 19:33:33Z spm $


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
