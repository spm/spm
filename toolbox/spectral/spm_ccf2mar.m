function [mar,pcond] = spm_ccf2mar(ccf,p)
% Converts cross covariance function to cross spectral density
% FORMAT [mar] = spm_ccf2mar(ccf,p)
%
% ccf  (N,m,m)   - cross covariance functions
% p              - AR(p) order
%
% mar.noise_cov  - (m,m)         covariance of innovations   
% mar.mean       - (p*m,m)       MAR coeficients (matrix format - positive)
% mar.lag        - lag(p).a(m,m) MAR coeficients (array format  - negative)
% mar.p          - order of a AR(p) model
%
% See also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_dcm_mtf.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_ccf2mar.m 6254 2014-11-04 18:24:21Z karl $


% MAR coeficients
%==========================================================================
N     = size(ccf,1);
m     = size(ccf,2);
n     = (N - 1)/2;
p     = min(p,n - 1);
ccf   = ccf((1:n) + n,:,:);

% create ccf matrices
%--------------------------------------------------------------------------
warning('off','MATLAB:toeplitz:DiagonalConflict')
A     = cell(m,m);
B     = cell(m,m);
for i = 1:m
    for j = 1:m
        A{i,j} = ccf((1:p) + 1,i,j);
        B{i,j} = toeplitz(ccf((1:p),i,j),ccf((1:p),j,i));
    end
end
warning('on','MATLAB:toeplitz:DiagonalConflict')



% least squares solution
%--------------------------------------------------------------------------
A    = spm_cat(A);
B    = spm_cat(B);
a    = full(B\A);

% convert mar (positive matrix) format to lag (negative array) format
%==========================================================================
mar.noise_cov  = squeeze(ccf(1,:,:)) - A'*a;
mar.lag        = spm_mar2lag(a); 
mar.a          = a;
mar.p          = p;
mar.d          = m;

% condition number of ccf matrix
%--------------------------------------------------------------------------
if nargout > 1
    pcond = cond(full(B));
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
