function [ccf] = spm_mar2ccf(mar,n)
% Get the cross covariance function from MAR coefficients or structure
% FORMAT [ccf] = spm_mar2ccf(mar,n)
%
% mar   - MAR coefficients or structure (see spm_mar.m)
% n     - number of time bins
%
% ccf   - (2*n + 1,i,j) cross covariance functions between I and J
%
% The mar coefficients are either specified in a cell array (as per
% spm_mar) or as a vector of (positive) coefficients as per spm_Q. The
% former are the negative values of the latter. If mar is a matrix of size
% d*p x d - it is assumed that the (positive) coefficients  run fast over 
% lag = p, as per the DCM routines.
%
% see also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________
% Copyright (C) 2014 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mar2ccf.m 6481 2015-06-16 17:01:47Z karl $


% Nyquist
%--------------------------------------------------------------------------
if nargin < 2, n = 128; end


% format coefficients into an array of negative coeficients (cf lag.a)
%--------------------------------------------------------------------------
if isvector(mar) && isnumeric(mar)
    mar = mar(:);
end
if isnumeric(mar)
    d  = size(mar,2);
    p  = size(mar,1)/d;
    for i = 1:d
        for j = 1:d
            for k = 1:p
                lag(k).a(i,j) = -mar((i - 1)*p + k,j);
            end
        end
    end
    mar = lag;
else
    d  = length(mar.lag(1).a);
    p  = length(mar.lag);
end

% covariance of innovations
%--------------------------------------------------------------------------
try
    c = mar.noise_cov;
catch
    c = eye(d,d);
end

% create AR representation and associated convolution kernels
%--------------------------------------------------------------------------
A     = cell(d,d);
B     = cell(d,d);
C     = cell(d,d);
for i = 1:d
    for j = 1:d
        a      = [0; mar.a((1:p) + (i - 1)*p,j)];
        A{i,j} = spdiags(ones(n,1)*a',-(0:p),n,n);
        C{i,j} = speye(n,n)*c(i,j);
    end
    B{i,i}     = speye(n,n);
end
A     = spm_cat(A);
B     = spm_cat(B);
C     = spm_cat(C);
K     = inv(B - A);

% compute cross-covariance matrices and reduces to an array of vectors
%--------------------------------------------------------------------------
CCF   = K*C*K';
ccf   = zeros(n,d,d);
for i = 1:d
    for j = 1:d        
        ccf(:,i,j) = full(CCF((1:n) + (i - 1)*n,ceil(n/2) + (j - 1)*n));
    end
end
