function [y,scalefactor] = spm_cond_units(y,n)
% Scales numeric arrays by a multiple of 10^n to avoid numerical overflow
% FORMAT [y,scalefector] = spm_cond_units(y,n)
%   y - y*scalefactor;
%   n - default 3
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_cond_units.m 3205 2009-06-16 10:15:00Z vladimir $
 
% default n = 3
%--------------------------------------------------------------------------
try, n; catch, n = 3; end

% rescale
%--------------------------------------------------------------------------
d           = spm_vec(y);
scalefactor = norm(d(~isnan(d)),1);
scalefactor = (10^n)^-round(log10(scalefactor)/n);
y           = spm_unvec(d*scalefactor,y);
