function [csd,dtf,coh,pha] = spm_mar2csd(mar,freqs,ns)
% Get spectral estimates from MAR model
% FORMAT [csd,dtf,coh,pha] = spm_mar2csd(mar,freqs,ns)
%
% mar   - MAR coefficients (see spm_mar.m)
% freqs - [Nf x 1] vector of frequencies to evaluate spectra at
% ns    - samples per second
%
% csd   - cross spectral density (assuming unit normal innovations)
% coh   - coherence
% pha   - phase
% dtf   - directed transfer function
%
% tthe man coefficients are either specified  as a cell array (as per
% spm_mar) or as a vector of (positive) coefficients as per spm_Q. The
% former are the negative values of the latter. If mar is a matrix of size
% d*p x d - it is assumed that the (positive) coefficients  run fast over 
% lag = p, as per the DCM routines.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Will Penny
% $Id: spm_mar_spectra.m 4096 2010-10-22 19:40:34Z karl $

% format coefficients
%--------------------------------------------------------------------------
if isvector(mar)
    mar = mar(:);
end
if ismatrix(mar)
    d  = size(mar,2);
    p  = size(mar,1)/d;
    for i = 1:d
        for j = 1:d
            for k = 1:p
                a(k).a(i,j) = -mar((i - 1)*p + k,j);
            end
        end
    end
    mar = a;
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
        af_tmp = af_tmp + mar(k).a*exp(-1i*w(ff)*k);
    end
    iaf_tmp     = inv(af_tmp);
    dtf(ff,:,:) = iaf_tmp;                            % transfer function
    csd(ff,:,:) = iaf_tmp*iaf_tmp';                   % and csd
end

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
