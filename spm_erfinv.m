function [u] = spm_erfinv(i)
% inverse of spm_erf (the error function)
% FORMAT [u] = spm_erfinv(i);
% i   - integral
% u   - threshold
%___________________________________________________________________________
% This routine is superceded by spm_invNcdf
%
%__________________________________________________________________________
% %W% %E%

fprintf('%cspm_erfinv is grandfathered, please replace with spm_invNcdf\n',7)

u = spm_invNcdf(1-i);

return

% Old way
%---------------------------------------------------------------------------


u = sqrt(2)*erfinv(1 - 2*i);
