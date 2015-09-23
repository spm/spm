function [csd,mtf,coh,pha] = spm_mar2csd(mar,freqs,ns)
% Get spectral estimates from MAR model
% FORMAT [csd,dtf,coh,pha] = spm_mar2csd(mar,freqs,ns)
%
% mar   - MAR coefficients or structure (see spm_mar.m)
% freqs - [Nf x 1] vector of frequencies to evaluate spectra at
% ns    - samples per second
%
% csd   - cross spectral density
% coh   - coherence
% pha   - phase
% mtf   - modulation transfer function
%
% The mar coefficients are either specified in a cell array (as per
% spm_mar) or as a vector of (positive) coefficients as per spm_Q. The
% former are the negative values of the latter. If mar is a matrix of size
% d*p x d - it is assumed that the (positive) coefficients  run fast over 
% lag = p, as per the DCM routines.
%
% see also:
%  spm_ccf2csd.m, spm_ccf2mar, spm_csd2ccf.m, spm_csd2mar.m, spm_mar2ccf.m,
%  spm_csd2coh.m, spm_Q.m, spm_mar.m and spm_mar_spectral.m
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mar2csd.m 6560 2015-09-23 13:50:43Z karl $


% Nyquist
%--------------------------------------------------------------------------
if nargin < 3, ns = 2*freqs(end); end

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
    clear mar
    mar.lag = lag;
else
    d  = length(mar.lag(1).a);
    p  = length(mar.lag);
end


% precision of innovation
%--------------------------------------------------------------------------
if isfield(mar,'noise_cov');
    C = mar.noise_cov;
else
    C = eye(d,d);
end


% frequencies
%--------------------------------------------------------------------------
Nf = length(freqs);
w  = 2*pi*freqs/ns;

% transfer function and complex cross spectral density
%--------------------------------------------------------------------------
for ff = 1:Nf,
    af_tmp = eye(d,d);
    for k = 1:p,
        af_tmp = af_tmp + mar.lag(k).a*exp(-1i*w(ff)*k);
    end
    iaf_tmp     = inv(af_tmp);
    mtf(ff,:,:) = iaf_tmp;                            % transfer function
    csd(ff,:,:) = iaf_tmp*C*iaf_tmp';                 % and csd
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
