function [i] = spm_erf(u)
% Returns the intergral from u to inf under the Unit Gaussian pdf
% FORMAT [i] = spm_erf(u);
% u   - threshold
% i   - integral
%___________________________________________________________________________
% This routine is superceded by spm_Ncdf
%
%__________________________________________________________________________
% %W% %E%

fprintf('%cspm_erf is grandfathered, please replace with spm_Ncdf\n',7)

i = 1 - spm_Ncdf(u);

return

% Old way
%---------------------------------------------------------------------------
i = (1 - erf(u/sqrt(2)))/2;
