function [ccf] = spm_mar2ccf(mar,n)
% gets the cross covariance function from MAR coefficients or structure
% FORMAT [ccf] = spm_mar2ccf(mar)
%
% mar   - MAR coefficients (see spm_mar.m)
% n     - number of time bins
%
% ccf   - (2*n + 1,i,j) cross covariance functions between I and J

%
% The mar coefficients are either specified as a cell array (as per
% spm_mar) or as a vector of (positive) coefficients as per spm_Q. The
% former are the negative values of the latter. If mar is a matrix of size
% d*p x d - it is assumed that the (positive) coefficients  run fast over 
% lag = p, as per the DCM routines.
%
% see also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2csd.m,
%  spm_csd2coh.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mar2csd.m 5889 2014-02-20 11:42:20Z karl $


% Nyquist
%--------------------------------------------------------------------------
if nargin < 2, n = 128; end

% format coefficients into an array of negative coeficients (cf lag.a)
%--------------------------------------------------------------------------
if isvector(mar)
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
    d  = length(mar(1).a);
    p  = length(mar);
end

% covariance innovations
%--------------------------------------------------------------------------
try
    C  = mar.noise_cov;
catch
    C  = eye(d,d);
end

% transfer function and complex cross spectral density
%--------------------------------------------------------------------------
for ff = 1:Nf,
    af_tmp = eye(d,d);
    for k = 1:p,
        af_tmp = af_tmp + mar(k).a*exp(-1i*w(ff)*k);
    end
    iaf_tmp     = inv(af_tmp);
    dtf(ff,:,:) = iaf_tmp;                            % transfer function
    csd(ff,:,:) = iaf_tmp*iaf_tmp';                   % and csd
end

% Normalise cross spectral density 
%--------------------------------------------------------------------------
csd = 2*csd/ns;

if nargout < 3, return, end

% Coherence and Phase
%--------------------------------------------------------------------------
for k = 1:d
    for j = 1:d
        rkj        = csd(:,k,j)./(sqrt(csd(:,k,k)).*sqrt(csd(:,j,j)));
        coh(:,k,j) = abs(rkj);
        pha(:,k,j) = atan(imag(rkj)./real(rkj));
    end
end
