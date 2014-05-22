function [Y,alpha] = spm_resample(X,alpha)
% Basic resample function (when no Signal Proc. Toolbox)
% FORMAT Y = spm_resample(X,alpha)
% IN:
%   - X: a nXm matrix of n time series
%   - alpha: the ration of input versus output sampling frequencies. If
%   alpha>1, rs(X,alpha) performs upsampling of the time series.
% OUT:
%   - Y: nX[alpha*m] matrix of resampled time series
%   - alpha: true alpha used (due to rational rounding)
% This function operates on rows of a signal matrix. This means it can be
% used on a block of channels.
%______________________________________________________________________________________
% Copyright (C) 2008-2014 Wellcome Trust Centre for Neuroimaging

% Jean Daunizeau

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
