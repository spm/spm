function [u] = spm_u(a,df,STAT)
% uncorrected critical height threshold at a specified significance level
% FORMAT [u] = spm_u(a,df,STAT)
% a     - critical probability - {alpha}
% df    - [df{interest} df{error}]
% STAT  - Statisical feild
%               'Z' - Gaussian feild
%               'T' - T - feild
%               'X' - Chi squared feild
%               'F' - F - feild
%
% u     - critical height {uncorrected}
%___________________________________________________________________________
% spm_u returns the uncorrected critical threshold at a specified significance
%
%___________________________________________________________________________
% Copyright (C) 2005 Wellcome Department of Imaging Neuroscience

% Karl Friston
% $Id: spm_u.m 112 2005-05-04 18:20:52Z john $


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
