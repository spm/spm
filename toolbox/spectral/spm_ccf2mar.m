function [mar,lag] = spm_ccf2mar(ccf,p)
% Converts cross covariance function to cross spectral density
% FORMAT [mar] = spm_ccf2mar(ccf,p)
%
% ccf  (N,m,m)          - cross covariance functions
% p                     - AR(p) order
%
% mar  (p*m,m)          - MAR coeficients (matrix format - positive)
% lag  lag(p).a(m,m)    - MAR coeficients (array format  - c.f., lag.a)
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ccf2mar.m 5837 2014-01-18 18:38:07Z karl $


% MAR coeficients
%==========================================================================
N     = size(ccf,1);
m     = size(ccf,2);
n     = (N - 1)/2;
ccf   = ccf((1:n) + n,:,:);

% create ccf matrices
%--------------------------------------------------------------------------
A     = cell(m,m);
B     = cell(m,m);
for i = 1:m
    for j = 1:m
        A{i,j} = ccf((1:p) + 1,i,j);
        B{i,j} = toeplitz(ccf((1:p),i,j));
    end
end

% least squares solution
%--------------------------------------------------------------------------
A    = spm_cat(A);
B    = spm_cat(B);
mar  = full(B\A);

% convert mar (positive matrix) format to lag (negative array) format
%==========================================================================
if nargout > 1
    lag = spm_mar2lag(mar);
end

function lag = spm_mar2lag(mar)
%--------------------------------------------------------------------------
[n m] = size(mar);
p     = n/m;
for i = 1:m
    for j = 1:m
        for k = 1:p
            lag(k).a(i,j) = -mar((i - 1)*p + k,j);
        end
    end
end
