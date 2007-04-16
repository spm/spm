function [p] = spm_wavspec (x,freqs,fs,show)
% Wavelet based spectrogram
% FORMAT [p] = spm_wavspec (x,freqs,fs,show)
% x         Data vector
% freqs     Frequencies to estimate power at
% fs        sample rate
% show      1 to plot real part of wavelet basis used (default = 0)
%___________________________________________________________________________
% Copyright (C) 2007 Wellcome Department of Imaging Neuroscience

% Will Penny 
% $Id$

x=x(:)';
freqs=freqs(:)';

if nargin < 4 | isempty(show)
    show=0;
end

Nf=length(freqs);
Nsamp=length(x);
M=spm_eeg_morlet(7,1000/fs,freqs);

if show
    dim=floor(sqrt(Nf));
    max_show=dim^2;
    figure
end

for i = 1:Nf,
    tmp = conv(x, M{i});
    
    if show & i < max_show+1
        subplot(dim,dim,i);
        t=[1:length(M{i})]/fs;
        plot(t,real(M{i}));
        title(sprintf('%1.2f',freqs(i)));
    end
        
    % time shift to remove delay
    tmp = tmp([1:Nsamp] + (length(M{i})-1)/2);
    p(i, :) = tmp.*conj(tmp);
end