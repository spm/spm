function [Y,alpha] = spm_timeseries_resample(X,alpha)
% Basic resample function (when no Signal Proc. Toolbox)
% FORMAT [Y,alpha] = spm_timeseries_resample(X,alpha)
% X      - (n x m) matrix of n time series of length m
% alpha  - the ratio of input versus output sampling frequencies.
%          If alpha>1, this performs upsampling of the time series.
%
% Y      - (n x [alpha*m]) matrix of resampled time series
% alpha  - true alpha used (due to rational rounding)
%
% This function operates on rows of a signal matrix.
%__________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau
% $Id: spm_timeseries_resample.m 6016 2014-05-23 17:34:06Z guillaume $

N0     = size(X,2);
N      = floor(N0*alpha);
alpha  = N/N0;
Y      = fftshift(fft(X,[],2),2);
sy     = size(Y,2);
middle = floor(sy./2)+1;
if alpha>1 % upsample
    N2 = floor((N-N0)./2);
    if N0/2 == floor(N0/2)
        Y(:,1) = []; % throw away non symmetric DFT coef
    end
    Y  = [zeros(size(Y,1),N2),Y,zeros(size(Y,1),N2)];
else % downsample
    N2 = floor(N./2);
    Y  = Y(:,middle-N2:middle+N2);
end
Y      = alpha*ifft(ifftshift(Y,2),[],2);
