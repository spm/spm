function M = spm_eeg_morlet(Rtf, ST, f, ff)
% Generate Morlet wavelets
% FORMAT M = spm_eeg_morlet(Rtf, ST, f, ff)
% 
% Rtf - 'wavelet factor', see [1]
% ST  - sample time [ms]
% f   - vector of frequencies [Hz]
% ff  - frequency to fix Gaussian envelope (sigma = Rtf/(2*pi*ff))
%       Default is ff = f, ie.e, a Morlet transform
%       NB: FWHM = sqrt(8*log(2))*sigma_t;
%
% M   - cell vector, where each element contains the filter for each
%       frequency in f
%__________________________________________________________________________
% 
% spm_eeg_morlet generates morlet wavelets for specified frequencies f with
% a specified ratio Rtf, see [1], for sample time ST (ms). One obtains the
% wavelet coefficients by convolution of a data vector with the kernels in
% M. See spm_eeg_tf how one obtains instantaneous power and phase estimates
% from the wavelet coefficients.
%
% [1] C. Tallon-Baudry, O. Bertrand, F. Peronnet and J. Pernier, 1998.
% Induced gamma-Band Activity during the Delay of a Visual Short-term
% memory Task in Humans. The Journal of Neuroscience (18): 4244-4254.
%__________________________________________________________________________
% Copyright (C) 2005-2017 Wellcome Trust Centre for Neuroimaging

% Stefan Kiebel
% $Id: spm_eeg_morlet.m 7122 2017-06-22 14:54:01Z guillaume $


if nargin < 4
    ff = f;
else
    ff = repmat(ff,1,numel(f));
end

M = cell(1,numel(f));

for i = 1:numel(f)
    
    % fixed or scale-dependent window
    %----------------------------------------------------------------------
    sigma_t = Rtf/(2*pi*ff(i));
    
    % this scaling factor is proportional to (Tallon-Baudry, 1998): 
    % (sigma_t*sqrt(pi))^(-1/2);
    %----------------------------------------------------------------------
    t = 0:ST*0.001:5*sigma_t;
    t = [-t(end:-1:2) t];
    M{i} = exp(-t.^2/(2*sigma_t^2)) .* exp(2 * 1i * pi * f(i) *t);    
    M{i} = M{i} ./ (sqrt(0.5*sum(real(M{i}).^2 + imag(M{i}).^2)));
    M{i} = M{i} - mean(M{i});
    
end
